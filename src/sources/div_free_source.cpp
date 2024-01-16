#include "div_free_source.hpp"
#include "helmholtz_projector.hpp"

namespace hephaestus
{

/*
Returns the Helmholtz (divergence-free) projection P(g) of a source g
defined by
(∇.P(g), q) = 0
where
P(g) = g - ∇Q         Q ∈ H1
P(g) = g + β∇p - ∇×M

via the weak form:
(g, ∇q) - (∇Q, ∇q) - <P(g).n, q> = 0
*/
DivFreeSource::DivFreeSource(const hephaestus::InputParameters & params)
  : src_coef_name(params.GetParam<std::string>("SourceName")),
    src_gf_name(params.GetParam<std::string>("SourceName")),
    hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
    h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")),
    potential_gf_name(
        params.GetOptionalParam<std::string>("PotentialName", std::string("_source_potential"))),
    solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
        "SolverOptions", hephaestus::InputParameters())),
    perform_helmholtz_projection(params.GetOptionalParam<bool>("HelmholtzProjection", true)),
    h_curl_mass(nullptr)
{
}

void
DivFreeSource::Init(hephaestus::GridFunctions & gridfunctions,
                    const hephaestus::FESpaces & fespaces,
                    hephaestus::BCMap & bc_map,
                    hephaestus::Coefficients & coefficients)
{
  H1FESpace_ = fespaces.Get(h1_fespace_name);
  if (H1FESpace_ == nullptr)
  {
    const std::string error_message = h1_fespace_name + " not found in fespaces when "
                                                        "creating DivFreeSource\n";
    mfem::mfem_error(error_message.c_str());
  }
  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name);
  if (HCurlFESpace_ == nullptr)
  {
    const std::string error_message = hcurl_fespace_name + " not found in fespaces when "
                                                           "creating DivFreeSource\n";
    mfem::mfem_error(error_message.c_str());
  }
  if (coefficients.vectors.Has(src_coef_name))
  {
    sourceVecCoef = coefficients.vectors.Get(src_coef_name);
  }
  else
  {
    MFEM_ABORT("SOURCE NOT FOUND");
  }

  // NB: Register must be false to avoid double-free.
  div_free_src_gf = std::make_unique<mfem::ParGridFunction>(HCurlFESpace_);
  gridfunctions.Register(src_gf_name, div_free_src_gf.get(), false);

  g = std::make_unique<mfem::ParGridFunction>(HCurlFESpace_);
  gridfunctions.Register("_user_source", g.get(), false);

  q_ = std::make_unique<mfem::ParGridFunction>(H1FESpace_);
  gridfunctions.Register(potential_gf_name, q_.get(), false);

  bc_map_ = &bc_map;
  gridfunctions_ = &gridfunctions;
  fespaces_ = &fespaces;

  buildHCurlMass();
}

void
DivFreeSource::buildHCurlMass()
{
  h_curl_mass = std::make_unique<mfem::ParBilinearForm>(HCurlFESpace_);
  h_curl_mass->AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();
}

void
DivFreeSource::Apply(mfem::ParLinearForm * lf)
{
  // Find an averaged representation of current density in H(curl)*
  g->ProjectCoefficient(*sourceVecCoef);
  mfem::ParLinearForm J(HCurlFESpace_);
  J.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*sourceVecCoef));
  J.Assemble();
  {
    mfem::HypreParMatrix M;
    mfem::Vector X, RHS;
    mfem::Array<int> ess_tdof_list;
    h_curl_mass->FormLinearSystem(ess_tdof_list, *g, J, M, X, RHS);

    DefaultGMRESSolver solver(solver_options, M);
    solver.Mult(RHS, X);

    h_curl_mass->RecoverFEMSolution(X, J, *g);
  }

  *div_free_src_gf = *g;

  if (perform_helmholtz_projection)
  {

    hephaestus::InputParameters projector_pars;
    projector_pars.SetParam("VectorGridFunctionName", src_gf_name);
    projector_pars.SetParam("ScalarGridFunctionName", potential_gf_name);
    projector_pars.SetParam("H1FESpaceName", h1_fespace_name);
    projector_pars.SetParam("HCurlFESpaceName", hcurl_fespace_name);

    hephaestus::BCMap bcs;

    hephaestus::HelmholtzProjector projector(projector_pars);
    projector.Project(*gridfunctions_, *fespaces_, bcs);
  }

  // Add divergence free source to target linear form
  h_curl_mass->Update();
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();

  h_curl_mass->AddMult(*div_free_src_gf, *lf, 1.0);
}

void
DivFreeSource::SubtractSource(mfem::ParGridFunction * gf)
{
  h_curl_mass->AddMult(*div_free_src_gf, *gf, -1.0);
}

} // namespace hephaestus

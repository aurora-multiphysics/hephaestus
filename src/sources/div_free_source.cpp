#include "div_free_source.hpp"
#include "helmholtz_projector.hpp"

namespace hephaestus {

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
DivFreeSource::DivFreeSource(const hephaestus::InputParameters &params)
    : src_coef_name(params.GetParam<std::string>("SourceName")),
      src_gf_name(params.GetParam<std::string>("SourceName")),
      hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
      h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")),
      potential_gf_name(params.GetOptionalParam<std::string>(
          "PotentialName", std::string("_source_potential"))),
      solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters())),
      perform_helmholtz_projection(
          params.GetOptionalParam<bool>("HelmholtzProjection", true)),
      a0(NULL), h_curl_mass(NULL), weakDiv_(NULL), grad(NULL), a0_solver(NULL) {
}

void DivFreeSource::Init(hephaestus::GridFunctions &gridfunctions,
                         const hephaestus::FESpaces &fespaces,
                         hephaestus::BCMap &bc_map,
                         hephaestus::Coefficients &coefficients) {
  H1FESpace_ = fespaces.Get(h1_fespace_name);
  if (H1FESpace_ == NULL) {
    const std::string error_message = h1_fespace_name +
                                      " not found in fespaces when "
                                      "creating DivFreeSource\n";
    mfem::mfem_error(error_message.c_str());
  }
  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name);
  if (HCurlFESpace_ == NULL) {
    const std::string error_message = hcurl_fespace_name +
                                      " not found in fespaces when "
                                      "creating DivFreeSource\n";
    mfem::mfem_error(error_message.c_str());
  }
  if (coefficients.vectors.Has(src_coef_name)) {
    sourceVecCoef = coefficients.vectors.Get(src_coef_name);
  } else {
    MFEM_ABORT("SOURCE NOT FOUND");
  }

  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  gridfunctions.Register(src_gf_name, div_free_src_gf, false);
  g = new mfem::ParGridFunction(HCurlFESpace_);
  gridfunctions.Register("_user_source", g, false);
  q_ = new mfem::ParGridFunction(H1FESpace_);
  gridfunctions.Register(potential_gf_name, q_, false);

  bc_map_ = &bc_map;
  gridfunctions_ = &gridfunctions;
  fespaces_ = &fespaces;

  this->buildH1Diffusion();
  this->buildHCurlMass();
  this->buildWeakDiv();
  this->buildGrad();

  // a0(p, p') = (β ∇ p, ∇ p')
  gDiv_ = new mfem::ParLinearForm(H1FESpace_);
  A0 = new mfem::HypreParMatrix;
  X0 = new mfem::Vector;
  B0 = new mfem::Vector;
}

void DivFreeSource::buildH1Diffusion() {
  if (a0 != NULL) {
    delete a0;
  }
  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a0->Assemble();
  a0->Finalize();
}

void DivFreeSource::buildHCurlMass() {
  if (h_curl_mass != NULL) {
    delete h_curl_mass;
  }
  h_curl_mass = new mfem::ParBilinearForm(HCurlFESpace_);
  h_curl_mass->AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();
}

void DivFreeSource::buildWeakDiv() {
  if (weakDiv_ != NULL) {
    delete weakDiv_;
  }
  weakDiv_ = new mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_);
  weakDiv_->AddDomainIntegrator(new mfem::VectorFEWeakDivergenceIntegrator);
  weakDiv_->Assemble();
  weakDiv_->Finalize();
}

void DivFreeSource::buildGrad() {
  if (grad != NULL) {
    delete grad;
  }
  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();
  grad->Finalize();
}

void DivFreeSource::Apply(mfem::ParLinearForm *lf) {
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

  if (perform_helmholtz_projection) {

    hephaestus::InputParameters projector_pars;
    projector_pars.SetParam("VectorGridFunctionName",
                             src_gf_name);
    projector_pars.SetParam("ScalarGridFunctionName",
                             potential_gf_name);
    projector_pars.SetParam("H1FESpaceName",
                             h1_fespace_name);
    projector_pars.SetParam("HCurlFESpaceName",
                             hcurl_fespace_name);

    hephaestus::HelmholtzProjector projector(projector_pars);
    projector.Project(*gridfunctions_, *fespaces_, *bc_map_);
  }

  // Add divergence free source to target linear form
  h_curl_mass->Update();
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();

  h_curl_mass->AddMult(*div_free_src_gf, *lf, 1.0);
}

void DivFreeSource::SubtractSource(mfem::ParGridFunction *gf) {
  h_curl_mass->AddMult(*div_free_src_gf, *gf, -1.0);
}

} // namespace hephaestus

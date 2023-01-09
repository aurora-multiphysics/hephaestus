#include "sources.hpp"

namespace hephaestus {

void Sources::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::DomainProperties &domain_properties) {
  for (const auto &[name, source] : GetMap()) {
    source->Init(variables, fespaces, domain_properties);
  }
}
void Sources::ApplySources(mfem::ParLinearForm *lf) {
  for (const auto &[name, source] : GetMap()) {
    source->ApplySource(lf);
  }
}
DivFreeVolumetricSource::DivFreeVolumetricSource(
    const hephaestus::InputParameters &params)
    : src_coef_name(params.GetParam<std::string>("SourceName")),
      src_gf_name(params.GetParam<std::string>("SourceName")),
      hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
      h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")),
      solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters())) {}

void DivFreeVolumetricSource::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::DomainProperties &domain_properties) {
  H1FESpace_ = fespaces.Get(h1_fespace_name);
  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name);

  if (domain_properties.vector_property_map.find(src_coef_name) !=
      domain_properties.vector_property_map.end()) {
    sourceVecCoef = domain_properties.vector_property_map[src_coef_name];
  } else {
    std::cout << "SOURCE NOT FOUND";
    exit;
  }

  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  variables.Register(src_gf_name, div_free_src_gf, false);
}

void DivFreeVolumetricSource::ApplySource(mfem::ParLinearForm *lf) {
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * 2;
  int geom = H1FESpace_->GetFE(0)->GetGeomType();
  const mfem::IntegrationRule *ir = &mfem::IntRules.Get(geom, irOrder);
  divFreeProj = new mfem::common::DivergenceFreeProjector(
      *H1FESpace_, *HCurlFESpace_, irOrder, NULL, NULL, NULL);

  /// Create a H(curl) mass matrix for integrating grid functions
  mfem::VectorFEMassIntegrator *h_curl_mass_integ =
      new mfem::VectorFEMassIntegrator;
  mfem::ParBilinearForm *h_curl_mass = new mfem::ParBilinearForm(HCurlFESpace_);

  h_curl_mass_integ->SetIntRule(ir);

  h_curl_mass->AddDomainIntegrator(h_curl_mass_integ);
  // assemble mass matrix
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();

  mfem::ParLinearForm J(HCurlFESpace_);
  J.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*sourceVecCoef));
  J.Assemble();

  mfem::ParGridFunction j(HCurlFESpace_);
  j.ProjectCoefficient(*sourceVecCoef);
  {
    mfem::HypreParMatrix M;
    mfem::Vector X, RHS;
    mfem::Array<int> ess_tdof_list;
    h_curl_mass->FormLinearSystem(ess_tdof_list, j, J, M, X, RHS);

    DefaultGMRESSolver solver(solver_options, M);
    solver.Mult(RHS, X);

    h_curl_mass->RecoverFEMSolution(X, J, j);
  }
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();
  *div_free_src_gf = 0.0;
  divFreeProj->Mult(j, *div_free_src_gf);

  // Compute the dual of div_free_src_gf
  h_curl_mass->AddMult(*div_free_src_gf, *lf);
}

} // namespace hephaestus

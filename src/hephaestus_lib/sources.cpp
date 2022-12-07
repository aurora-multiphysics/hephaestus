#include "sources.hpp"

namespace hephaestus {

DivFreeVolumetricSource::DivFreeVolumetricSource(
    const hephaestus::InputParameters &params)
    : src_gf_name(params.GetParam<std::string>("SourceName")),
      hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
      h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")) {}

void DivFreeVolumetricSource::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces) {
  H1FESpace_ = fespaces.Get(h1_fespace_name);
  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name);

  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  // _variables.Register(src_gf_name, div_free_src_gf, false);

  /// get int rule (approach followed my MFEM Tesla Miniapp)
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * 2;
  int geom = H1FESpace_->GetFE(0)->GetGeomType();
  const mfem::IntegrationRule *ir = &mfem::IntRules.Get(geom, irOrder);
  divFreeProj = new mfem::common::DivergenceFreeProjector(
      *H1FESpace_, *HCurlFESpace_, irOrder, NULL, NULL, NULL);

  /// Create a H(curl) mass matrix for integrating grid functions
  mfem::BilinearFormIntegrator *h_curl_mass_integ =
      new mfem::VectorFEMassIntegrator;
  h_curl_mass_integ->SetIntRule(ir);

  h_curl_mass = new mfem::ParBilinearForm(HCurlFESpace_);
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

    mfem::HypreBoomerAMG amg(M);
    amg.SetPrintLevel(0);
    mfem::HypreGMRES gmres(M);
    gmres.SetTol(1e-12);
    gmres.SetMaxIter(200);
    gmres.SetPrintLevel(0);
    gmres.SetPreconditioner(amg);
    gmres.Mult(RHS, X);

    h_curl_mass->RecoverFEMSolution(X, J, j);
  }
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();

  *div_free_src_gf = 0.0;
  divFreeProj->Mult(j, *div_free_src_gf);
}

void DivFreeVolumetricSource::ApplySource(mfem::LinearForm *lf) {
  // Compute the dual of div_free_src_gf
  h_curl_mass->AddMult(*div_free_src_gf, *lf);
}

} // namespace hephaestus

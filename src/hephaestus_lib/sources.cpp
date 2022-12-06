#include "sources.hpp"

namespace hephaestus {

DivFreeVolumetricSource::DivFreeVolumetricSource(
    mfem::VectorFunctionCoefficient *sourceVecCoef) {

  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  // _variables.Register("source", div_free_src_gf, false);
  // int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() +
  //               2 * H1FESpace_->GetOrder();
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * 2;
  divFreeProj = new mfem::common::DivergenceFreeProjector(
      *H1FESpace_, *HCurlFESpace_, irOrder, NULL, NULL, NULL);
  hCurlMass = new mfem::ParBilinearForm(HCurlFESpace_);
  hCurlMass->AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  hCurlMass->Assemble();

  /// get int rule (approach followed my MFEM Tesla Miniapp)
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * 2;
  int geom = H1FESpace_->GetFE(0)->GetGeomType();
  const mfem::IntegrationRule *ir = &mfem::IntRules.Get(geom, irOrder);

  /// Create a H(curl) mass matrix for integrating grid functions
  mfem::BilinearFormIntegrator *h_curl_mass_integ =
      new mfem::VectorFEMassIntegrator;
  h_curl_mass_integ->SetIntRule(ir);
  mfem::ParBilinearForm h_curl_mass(HCurlFESpace_);
  h_curl_mass.AddDomainIntegrator(h_curl_mass_integ);
  // assemble mass matrix
  h_curl_mass.Assemble();
  h_curl_mass.Finalize();

  mfem::ParLinearForm J(HCurlFESpace_);
  J.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*sourceVecCoef));
  J.Assemble();

  mfem::ParGridFunction j(HCurlFESpace_);
  j.ProjectCoefficient(*sourceVecCoef);
  {
    mfem::HypreParMatrix M;
    mfem::Vector X, RHS;
    mfem::Array<int> ess_tdof_list;
    h_curl_mass.FormLinearSystem(ess_tdof_list, j, J, M, X, RHS);

    // if (rank == 0)
    //   std::cout << "solving for J in H(curl)\n";

    mfem::HypreBoomerAMG amg(M);
    amg.SetPrintLevel(0);
    mfem::HypreGMRES gmres(M);
    gmres.SetTol(1e-12);
    gmres.SetMaxIter(200);
    gmres.SetPrintLevel(0);
    gmres.SetPreconditioner(amg);
    gmres.Mult(RHS, X);

    h_curl_mass.RecoverFEMSolution(X, J, j);
  }
  h_curl_mass.Assemble();
  h_curl_mass.Finalize();

  *div_free_src_gf = 0.0;
  divFreeProj->Mult(j, *div_free_src_gf);

  // src_gf->ProjectCoefficient(*sourceVecCoef);
  // Compute the discretely divergence-free portion of src_gf
  // divFreeProj->Mult(*src_gf, *div_free_src_gf);
  // Compute the dual of div_free_src_gf
  h_curl_mass.AddMult(*div_free_src_gf, *b1);
}

} // namespace hephaestus

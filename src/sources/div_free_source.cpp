#include "div_free_source.hpp"

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

(-Q u, grad v)
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
      a0(NULL), h_curl_mass(NULL), weakDiv_(NULL), grad(NULL), a0_solver(NULL) {
}

void DivFreeSource::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
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
  if (domain_properties.vector_property_map.Has(src_coef_name)) {
    sourceVecCoef = domain_properties.vector_property_map.Get(src_coef_name);
  } else {
    MFEM_ABORT("SOURCE NOT FOUND");
  }

  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  variables.Register(src_gf_name, div_free_src_gf, false);
  q_ = new mfem::ParGridFunction(H1FESpace_);
  variables.Register(potential_gf_name, q_, false);

  _bc_map = &bc_map;

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
  grad->Assemble();
  grad->Finalize();
}

void DivFreeSource::Apply(mfem::ParLinearForm *lf) {
  // (g, ∇q) - (∇Q, ∇q) - <P(g).n, q> = 0
  mfem::ParGridFunction g(HCurlFESpace_);
  g.ProjectCoefficient(*sourceVecCoef);

  mfem::ParLinearForm J(HCurlFESpace_);
  J.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*sourceVecCoef));
  J.Assemble();

  {
    mfem::HypreParMatrix M;
    mfem::Vector X, RHS;
    mfem::Array<int> ess_tdof_list;
    h_curl_mass->FormLinearSystem(ess_tdof_list, g, J, M, X, RHS);

    DefaultGMRESSolver solver(solver_options, M);
    solver.Mult(RHS, X);

    h_curl_mass->RecoverFEMSolution(X, J, g);
  }

  // begin div free proj
  *q_ = 0.0;
  int myid = H1FESpace_->GetMyRank();
  mfem::Array<int> ess_bdr_;
  mfem::Array<int> ess_bdr_tdofs_;
  mfem::ParGridFunction Phi_gf(H1FESpace_);
  // <P(g).n, q>
  _bc_map->applyEssentialBCs(potential_gf_name, ess_bdr_tdofs_, Phi_gf,
                             (H1FESpace_->GetParMesh()));
  _bc_map->applyIntegratedBCs(potential_gf_name, *gDiv_,
                              (H1FESpace_->GetParMesh()));
  gDiv_->Assemble();
  // Compute the divergence of g
  // // (g, ∇q)
  weakDiv_->AddMult(g, *gDiv_, -1.0);

  // Apply essential BC and form linear system
  ess_bdr_.SetSize(H1FESpace_->GetParMesh()->bdr_attributes.Max());
  ess_bdr_ = 0;
  ess_bdr_tdofs_.SetSize((myid == 0) ? 1 : 0);
  if (myid == 0) {
    ess_bdr_tdofs_[0] = 0;
  }

  // (g, ∇q) - (∇Q, ∇q) - <P(g).n, q> = 0
  // (∇Q, ∇q) = (g, ∇q) - <P(g).n, q>
  a0->FormLinearSystem(ess_bdr_tdofs_, *q_, *gDiv_, *A0, *X0, *B0);

  // Solve the linear system for Psi
  mfem::HypreBoomerAMG *amg_ = new mfem::HypreBoomerAMG(*A0);
  amg_->SetPrintLevel(0);
  mfem::HyprePCG *pcg_ = new mfem::HyprePCG(*A0);
  pcg_->SetTol(1e-14);
  pcg_->SetMaxIter(200);
  pcg_->SetPrintLevel(0);
  pcg_->SetPreconditioner(*amg_);
  pcg_->Mult(*B0, *X0);
  delete amg_;
  delete pcg_;

  a0->RecoverFEMSolution(*X0, *gDiv_, *q_);
  // Compute the irrotational component of g
  // P(g) = g - ∇Q
  grad->Mult(*q_, *div_free_src_gf);
  *div_free_src_gf -= g;
  *div_free_src_gf *= -1.0;
  // end div free proj
  // Compute the dual of div_free_src_gf
  h_curl_mass->Update();
  h_curl_mass->Assemble();
  h_curl_mass->Finalize();

  h_curl_mass->AddMult(*div_free_src_gf, *lf, 1.0);
}

void DivFreeSource::SubtractSource(mfem::ParGridFunction *gf) {
  h_curl_mass->AddMult(*div_free_src_gf, *gf, -1.0);
}

} // namespace hephaestus

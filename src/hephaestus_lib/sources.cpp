#include "sources.hpp"

namespace hephaestus {

void Sources::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  for (const auto &[name, source] : GetMap()) {
    source->Init(variables, fespaces, bc_map, domain_properties);
  }
}
void Sources::ApplySources(mfem::ParLinearForm *lf) {
  for (const auto &[name, source] : GetMap()) {
    source->ApplySource(lf);
  }
}

// Create a scalar potential source that can add terms of the form
// J0_{n+1} ∈ H(div) source field, where J0 = -β∇p and β is a conductivity
// coefficient.
ScalarPotentialSource::ScalarPotentialSource(
    const hephaestus::InputParameters &params)
    : src_gf_name(params.GetParam<std::string>("SourceName")),
      potential_gf_name(params.GetParam<std::string>("PotentialName")),
      hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
      h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")),
      beta_coef_name(params.GetParam<std::string>("ConductivityCoefName")),
      solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters())),
      grad(NULL), m1(NULL), a0_solver(NULL) {}

void ScalarPotentialSource::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  H1FESpace_ = fespaces.Get(h1_fespace_name);
  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name);

  p_ = new mfem::ParGridFunction(H1FESpace_);
  variables.Register(potential_gf_name, p_, false);
  p_ = variables.Get(potential_gf_name);

  _bc_map = &bc_map;

  betaCoef = domain_properties.scalar_property_map[beta_coef_name];

  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*betaCoef));
  a0->Assemble();

  this->buildGrad();
  this->buildM1(betaCoef);
  // a0(p, p') = (β ∇ p, ∇ p')
  b0 = new mfem::ParLinearForm(H1FESpace_);
  A0 = new mfem::HypreParMatrix;
  X0 = new mfem::Vector;
  B0 = new mfem::Vector;

  grad_p_ = new mfem::ParGridFunction(HCurlFESpace_);
  variables.Register(src_gf_name, grad_p_, false);
}

void ScalarPotentialSource::buildM1(mfem::Coefficient *Sigma) {
  if (m1 != NULL) {
    delete m1;
  }

  m1 = new mfem::ParBilinearForm(HCurlFESpace_);
  m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  m1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void ScalarPotentialSource::buildGrad() {
  if (grad != NULL) {
    delete grad;
  }

  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void ScalarPotentialSource::ApplySource(mfem::ParLinearForm *lf) {
  // -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
  // a0(p_{n+1}, p') = b0(p')
  // a0(p, p') = (β ∇ p, ∇ p')
  // b0(p') = <n.s0, p'>
  mfem::ParGridFunction Phi_gf(H1FESpace_);
  mfem::Array<int> poisson_ess_tdof_list;
  Phi_gf = 0.0;
  *b0 = 0.0;
  _bc_map->applyEssentialBCs(potential_gf_name, poisson_ess_tdof_list, Phi_gf,
                             (H1FESpace_->GetParMesh()));
  _bc_map->applyIntegratedBCs(potential_gf_name, *b0,
                              (H1FESpace_->GetParMesh()));
  b0->Assemble();

  a0->FormLinearSystem(poisson_ess_tdof_list, Phi_gf, *b0, *A0, *X0, *B0);

  if (a0_solver == NULL) {
    hephaestus::InputParameters h1_solver_options;
    h1_solver_options.SetParam("Tolerance", float(1.0e-9));
    h1_solver_options.SetParam("MaxIter", (unsigned int)1000);
    h1_solver_options.SetParam("PrintLevel", -1);
    a0_solver = new hephaestus::DefaultH1PCGSolver(h1_solver_options, *A0);
  }
  // Solve
  a0_solver->Mult(*B0, *X0);

  // "undo" the static condensation saving result in grid function dP
  a0->RecoverFEMSolution(*X0, *b0, *p_);

  grad->Mult(*p_, *grad_p_);
  m1->AddMult(*grad_p_, *lf, 1.0);
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
    hephaestus::BCMap &bc_map,
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

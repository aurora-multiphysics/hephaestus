#include "scalar_potential_source.hpp"

namespace hephaestus {

// Create a scalar potential source that can add terms of the form
// J0_{n+1} ∈ H(div) source field, where J0 = -β∇p and β is a conductivity
// coefficient.
ScalarPotentialSource::ScalarPotentialSource(
    const hephaestus::InputParameters &params)
    : efield_gf_name(params.GetParam<std::string>("EFieldName")),
      potential_gf_name(params.GetParam<std::string>("PotentialName")),
      hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
      h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")),
      beta_coef_name(params.GetParam<std::string>("ConductivityCoefName")),
      solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters())),
      grad(NULL), m1(NULL), a0_solver(NULL) {}

void ScalarPotentialSource::Init(hephaestus::GridFunctions &gridfunctions,
                                 const hephaestus::FESpaces &fespaces,
                                 hephaestus::BCMap &bc_map,
                                 hephaestus::Coefficients &coefficients) {
  H1FESpace_ = fespaces.Get(h1_fespace_name);
  if (H1FESpace_ == NULL) {
    const std::string error_message = h1_fespace_name +
                                      " not found in fespaces when "
                                      "creating ScalarPotentialSource\n";
    mfem::mfem_error(error_message.c_str());
  }

  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name);
  if (HCurlFESpace_ == NULL) {
    const std::string error_message = hcurl_fespace_name +
                                      " not found in fespaces when "
                                      "creating ScalarPotentialSource\n";
    mfem::mfem_error(error_message.c_str());
  }

  p_ = gridfunctions.Get(potential_gf_name);
  if (p_ == NULL) {
    std::cout << potential_gf_name + " not found in gridfunctions when "
                                     "creating ScalarPotentialSource. "
                                     "Creating new ParGridFunction\n";
    p_ = new mfem::ParGridFunction(H1FESpace_);
    gridfunctions.Register(potential_gf_name, p_, false);
  }

  grad_p_ = gridfunctions.Get(efield_gf_name);
  if (grad_p_ == NULL) {
    grad_p_ = new mfem::ParGridFunction(HCurlFESpace_);
    gridfunctions.Register(efield_gf_name, grad_p_, false);
  }

  _bc_map = &bc_map;

  betaCoef = coefficients.scalars.Get(beta_coef_name);

  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*betaCoef));
  a0->Assemble();

  buildGrad();
  buildM1(betaCoef);
  // a0(p, p') = (β ∇ p, ∇ p')
  b0 = new mfem::ParLinearForm(H1FESpace_);
  A0 = new mfem::HypreParMatrix;
  X0 = new mfem::Vector;
  B0 = new mfem::Vector;
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

void ScalarPotentialSource::Apply(mfem::ParLinearForm *lf) {
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

  a0->Update();
  a0->Assemble();
  a0->FormLinearSystem(poisson_ess_tdof_list, Phi_gf, *b0, *A0, *X0, *B0);

  if (a0_solver == NULL) {
    a0_solver = new hephaestus::DefaultH1PCGSolver(solver_options, *A0);
  }
  // Solve
  a0_solver->Mult(*B0, *X0);

  // "undo" the static condensation saving result in grid function dP
  a0->RecoverFEMSolution(*X0, *b0, *p_);

  grad->Mult(*p_, *grad_p_);

  m1->Update();
  m1->Assemble();
  m1->AddMult(*grad_p_, *lf, 1.0);
}

void ScalarPotentialSource::SubtractSource(mfem::ParGridFunction *gf) {
  grad->AddMult(*p_, *gf, -1.0);
}

} // namespace hephaestus

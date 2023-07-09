// Solves the equations
// ∇⋅s0 = 0
// ∇×(α∇×u) + βdu/dt = s0

// where
// s0 ∈ H(div) source field
// u ∈ H(curl)
// p ∈ H1

// Dirichlet boundaries constrain du/dt
// Integrated boundaries constrain (α∇×u) × n

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (α∇×u, ∇×u') + (βdu/dt, u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// du/dt_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
// - (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
// using
// u_{n+1} = u_{n} + dt du/dt_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(du/dt_{n+1}, u') = b1(u')

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>
#include "hcurl_formulation.hpp"

namespace hephaestus {

HCurlFormulation::HCurlFormulation() : TimeDomainFormulation() {
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
  h_curl_var_name = std::string("h_curl_var");
}

void HCurlFormulation::ConstructEquationSystem() {
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("HCurlVarName", h_curl_var_name);
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  this->GetProblem()->td_equation_system =
      std::make_unique<hephaestus::CurlCurlEquationSystem>(weak_form_params);
}

void HCurlFormulation::ConstructOperator() {
  this->problem->td_operator = std::make_unique<hephaestus::HCurlOperator>(
      *(this->problem->pmesh), this->problem->fespaces,
      this->problem->gridfunctions, this->problem->bc_map,
      this->problem->coefficients, this->problem->sources,
      this->problem->solver_options);
  this->problem->td_operator->SetEquationSystem(
      this->problem->td_equation_system.get());
  this->problem->td_operator->Setgridfunctions();
};

void HCurlFormulation::RegisterGridFunctions() {
  int &myid = this->GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::FESpaces &fespaces = this->GetProblem()->fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(h_curl_var_name)) {
    if (myid == 0) {
      MFEM_WARNING(h_curl_var_name << " not found in gridfunctions: building "
                                      "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(h_curl_var_name, std::string("_HCurlFESpace"));
  };
  // Register time derivatives
  TimeDomainProblemBuilder::RegisterGridFunctions();
};

CurlCurlEquationSystem::CurlCurlEquationSystem(
    const hephaestus::InputParameters &params)
    : TimeDependentEquationSystem(params),
      h_curl_var_name(params.GetParam<std::string>("HCurlVarName")),
      alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
      beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
      dtalpha_coef_name(std::string("dt_") + alpha_coef_name) {}

void CurlCurlEquationSystem::Init(hephaestus::GridFunctions &gridfunctions,
                                  const hephaestus::FESpaces &fespaces,
                                  hephaestus::BCMap &bc_map,
                                  hephaestus::Coefficients &coefficients) {
  coefficients.scalars.Register(
      dtalpha_coef_name,
      new mfem::TransformedCoefficient(
          &dtCoef, coefficients.scalars.Get(alpha_coef_name), prodFunc),
      true);
  TimeDependentEquationSystem::Init(gridfunctions, fespaces, bc_map,
                                    coefficients);
}

void CurlCurlEquationSystem::addKernels() {
  addVariableNameIfMissing(h_curl_var_name);
  std::string dh_curl_var_dt = GetTimeDerivativeName(h_curl_var_name);

  // (α∇×u_{n}, ∇×u')
  hephaestus::InputParameters weakCurlCurlParams;
  weakCurlCurlParams.SetParam("CoupledVariableName", h_curl_var_name);
  weakCurlCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(dh_curl_var_dt,
            new hephaestus::WeakCurlCurlKernel(weakCurlCurlParams));

  // (αdt∇×du/dt_{n+1}, ∇×u')
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(dh_curl_var_dt, new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βdu/dt_{n+1}, u')
  hephaestus::InputParameters vectorFEMassParams;
  vectorFEMassParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(dh_curl_var_dt,
            new hephaestus::VectorFEMassKernel(vectorFEMassParams));
}

void HCurlFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has(alpha_coef_name)) {
    MFEM_ABORT(alpha_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has(beta_coef_name)) {
    MFEM_ABORT(beta_coef_name + " coefficient not found.");
  }
}

HCurlOperator::HCurlOperator(mfem::ParMesh &pmesh,
                             hephaestus::FESpaces &fespaces,
                             hephaestus::GridFunctions &gridfunctions,
                             hephaestus::BCMap &bc_map,
                             hephaestus::Coefficients &coefficients,
                             hephaestus::Sources &sources,
                             hephaestus::InputParameters &solver_options)
    : TimeDomainEquationSystemOperator(pmesh, fespaces, gridfunctions, bc_map,
                                       coefficients, sources, solver_options) {}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing p, u and v.

Unknowns
s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
du/dt_{n+1} ∈ H(curl)
p_{n+1} ∈ H1

Fully discretised equations
-(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
(α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
- (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
using
u_{n+1} = u_{n} + dt du/dt_{n+1}
*/
void HCurlOperator::ImplicitSolve(const double dt, const mfem::Vector &X,
                                  mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    local_trial_vars.at(ind)->MakeRef(local_trial_vars.at(ind)->ParFESpace(),
                                      dX_dt, true_offsets[ind]);
  }
  _coefficients.SetTime(this->GetTime());
  _equation_system->setTimeStep(dt);
  _equation_system->updateEquationSystem(_bc_map, _sources);

  _equation_system->FormLinearSystem(blockA, trueX, trueRhs);

  if (a1_solver != NULL) {
    delete a1_solver;
  }
  a1_solver = new hephaestus::DefaultHCurlPCGSolver(
      _solver_options, *blockA.As<mfem::HypreParMatrix>(),
      _equation_system->test_pfespaces.at(0));
  a1_solver->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _gridfunctions);
}

} // namespace hephaestus

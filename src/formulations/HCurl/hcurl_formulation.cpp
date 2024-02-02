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

#include <utility>

namespace hephaestus
{

HCurlFormulation::HCurlFormulation(std::string alpha_coef_name,
                                   std::string beta_coef_name,
                                   std::string h_curl_var_name)
  : _alpha_coef_name(std::move(alpha_coef_name)),
    _beta_coef_name(std::move(beta_coef_name)),
    _h_curl_var_name(std::move(h_curl_var_name))
{
}

void
HCurlFormulation::ConstructEquationSystem()
{
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("HCurlVarName", _h_curl_var_name);
  weak_form_params.SetParam("AlphaCoefName", _alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", _beta_coef_name);
  GetProblem()->_td_equation_system =
      std::make_unique<hephaestus::CurlCurlEquationSystem>(weak_form_params);
}

void
HCurlFormulation::ConstructOperator()
{
  // Default solver options for the HCurl formulation
  hephaestus::InputParameters & solver_options = GetProblem()->_solver_options;
  solver_options.SetParamIfEmpty("LinearSolver", std::string("FGMRES"));
  solver_options.SetParamIfEmpty("LinearPreconditioner", std::string("AMS"));
  solver_options.SetParamIfEmpty("Tolerance", float(1.0e-13));
  solver_options.SetParamIfEmpty("AbsTolerance", float(1.0e-20));
  solver_options.SetParamIfEmpty("MaxIter", (unsigned int)500);
  solver_options.SetParamIfEmpty("PrintLevel", 1);
  _problem->_solvers.SetSolverOptions(solver_options);

  _problem->_td_operator = std::make_unique<hephaestus::HCurlOperator>(*(_problem->_pmesh),
                                                                       _problem->_fespaces,
                                                                       _problem->_gridfunctions,
                                                                       _problem->_bc_map,
                                                                       _problem->_coefficients,
                                                                       _problem->_sources,
                                                                       _problem->_solver_options,
                                                                       _problem->_solvers);
  _problem->_td_operator->SetEquationSystem(_problem->_td_equation_system.get());
  _problem->_td_operator->SetGridFunctions();
};

void
HCurlFormulation::RegisterGridFunctions()
{
  int & myid = GetProblem()->_myid;
  hephaestus::GridFunctions & gridfunctions = GetProblem()->_gridfunctions;
  hephaestus::FESpaces & fespaces = GetProblem()->_fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_h_curl_var_name))
  {
    if (myid == 0)
    {
      MFEM_WARNING(_h_curl_var_name << " not found in gridfunctions: building "
                                       "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(_h_curl_var_name, std::string("_HCurlFESpace"));
  };
  // Register time derivatives
  TimeDomainProblemBuilder::RegisterGridFunctions();
};

CurlCurlEquationSystem::CurlCurlEquationSystem(const hephaestus::InputParameters & params)
  : TimeDependentEquationSystem(params),
    _h_curl_var_name(params.GetParam<std::string>("HCurlVarName")),
    _alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
    _beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
    _dtalpha_coef_name(std::string("dt_") + _alpha_coef_name)
{
}

void
CurlCurlEquationSystem::Init(hephaestus::GridFunctions & gridfunctions,
                             const hephaestus::FESpaces & fespaces,
                             hephaestus::BCMap & bc_map,
                             hephaestus::Coefficients & coefficients)
{
  coefficients._scalars.Register(
      _dtalpha_coef_name,
      new mfem::TransformedCoefficient(
          &_dt_coef, coefficients._scalars.Get(_alpha_coef_name), prodFunc),
      true);
  TimeDependentEquationSystem::Init(gridfunctions, fespaces, bc_map, coefficients);
}

void
CurlCurlEquationSystem::AddKernels()
{
  AddVariableNameIfMissing(_h_curl_var_name);
  std::string dh_curl_var_dt = GetTimeDerivativeName(_h_curl_var_name);

  // (α∇×u_{n}, ∇×u')
  hephaestus::InputParameters weak_curl_curl_params;
  weak_curl_curl_params.SetParam("CoupledVariableName", _h_curl_var_name);
  weak_curl_curl_params.SetParam("CoefficientName", _alpha_coef_name);
  AddKernel(dh_curl_var_dt,
            std::make_shared<hephaestus::WeakCurlCurlKernel>(weak_curl_curl_params));

  // (αdt∇×du/dt_{n+1}, ∇×u')
  hephaestus::InputParameters curl_curl_params;
  curl_curl_params.SetParam("CoefficientName", _dtalpha_coef_name);
  AddKernel(dh_curl_var_dt, std::make_shared<hephaestus::CurlCurlKernel>(curl_curl_params));

  // (βdu/dt_{n+1}, u')
  hephaestus::InputParameters vector_fe_mass_params;
  vector_fe_mass_params.SetParam("CoefficientName", _beta_coef_name);
  AddKernel(dh_curl_var_dt,
            std::make_shared<hephaestus::VectorFEMassKernel>(vector_fe_mass_params));
}

void
HCurlFormulation::RegisterCoefficients()
{
  hephaestus::Coefficients & coefficients = GetProblem()->_coefficients;
  if (!coefficients._scalars.Has(_alpha_coef_name))
  {
    MFEM_ABORT(_alpha_coef_name + " coefficient not found.");
  }
  if (!coefficients._scalars.Has(_beta_coef_name))
  {
    MFEM_ABORT(_beta_coef_name + " coefficient not found.");
  }
}

HCurlOperator::HCurlOperator(mfem::ParMesh & pmesh,
                             hephaestus::FESpaces & fespaces,
                             hephaestus::GridFunctions & gridfunctions,
                             hephaestus::BCMap & bc_map,
                             hephaestus::Coefficients & coefficients,
                             hephaestus::Sources & sources,
                             hephaestus::InputParameters & solver_options,
                             hephaestus::ProblemSolvers & solvers)
  : TimeDomainEquationSystemOperator(
        pmesh, fespaces, gridfunctions, bc_map, coefficients, sources, solver_options),
    _solvers(&solvers)
{
}

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
void
HCurlOperator::ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt)
{
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _local_test_vars.at(ind)->MakeRef(
        _local_test_vars.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    _local_trial_vars.at(ind)->MakeRef(
        _local_trial_vars.at(ind)->ParFESpace(), dX_dt, _true_offsets[ind]);
  }
  _coefficients.SetTime(GetTime());
  _equation_system->SetTimeStep(dt);
  _equation_system->UpdateEquationSystem(_bc_map, _sources);

  _equation_system->FormLinearSystem(_block_a, _true_x, _true_rhs);

  _solvers->SetEdgeFESpace(_equation_system->_test_pfespaces.at(0));
  _solvers->SetLinearPreconditioner(*_block_a.As<mfem::HypreParMatrix>());
  _solvers->SetLinearSolver(*_block_a.As<mfem::HypreParMatrix>());
  _solvers->_linear_solver->Mult(_true_rhs, _true_x);

  _equation_system->RecoverFEMSolution(_true_x, _gridfunctions);
}

} // namespace hephaestus

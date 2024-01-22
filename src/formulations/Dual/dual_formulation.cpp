// Solves the equations
// ∇⋅s0 = 0
// ∇×(αv) - βu = s0
// dv/dt = -∇×u

// where
// s0 ∈ H(div) source field
// v ∈ H(div)
// u ∈ H(curl)
// p ∈ H1

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (αv, ∇×u') - (βu, u') - (s0, u') - <(αv) × n, u'> = 0
// (dv/dt, v') + (∇×u, v') = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// dv/dt_{n+1} ∈ H(div)
// u_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (αv_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
// <(αv) × n, u'> = 0
// (dv/dt_{n+1}, v') + (∇×u_{n+1}, v') = 0
// using
// v_{n+1} = v_{n} + dt dv/dt_{n+1} = v_{n} - dt ∇×u_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(u_{n+1}, u') = b1(u')
// dv/dt_{n+1} = -∇×u

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') + (αv_{n}, ∇×u') + <(αdt∇×u_{n+1}) × n, u'>

#include "dual_formulation.hpp"

#include <utility>

namespace hephaestus
{

DualFormulation::DualFormulation(std::string alpha_coef_name,
                                 std::string beta_coef_name,
                                 std::string h_curl_var_name,
                                 std::string h_div_var_name)
  : _alpha_coef_name(std::move(alpha_coef_name)),
    _beta_coef_name(std::move(beta_coef_name)),
    _h_curl_var_name(std::move(h_curl_var_name)),
    _h_div_var_name(std::move(h_div_var_name))
{
}

void
DualFormulation::ConstructEquationSystem()
{
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("HCurlVarName", _h_curl_var_name);
  weak_form_params.SetParam("HDivVarName", _h_div_var_name);
  weak_form_params.SetParam("AlphaCoefName", _alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", _beta_coef_name);
  GetProblem()->_td_equation_system =
      std::make_unique<hephaestus::WeakCurlEquationSystem>(weak_form_params);
}

void
DualFormulation::ConstructOperator()
{
  _problem->_td_operator = std::make_unique<hephaestus::DualOperator>(*(_problem->_pmesh),
                                                                      _problem->_fespaces,
                                                                      _problem->_gridfunctions,
                                                                      _problem->_bc_map,
                                                                      _problem->_coefficients,
                                                                      _problem->_sources,
                                                                      _problem->_solver_options);
  _problem->_td_operator->SetEquationSystem(_problem->_td_equation_system.get());
  _problem->_td_operator->SetGridFunctions();
};

void
DualFormulation::RegisterGridFunctions()
{
  int & myid = GetProblem()->_myid;
  hephaestus::GridFunctions & gridfunctions = GetProblem()->_gridfunctions;

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
  }

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_h_div_var_name))
  {
    if (myid == 0)
    {
      MFEM_WARNING(_h_div_var_name << " not found in gridfunctions: building "
                                      "gridfunction from defaults");
    }
    AddFESpace(std::string("_HDivFESpace"), std::string("RT_3D_P3"));
    AddGridFunction(_h_div_var_name, std::string("_HDivFESpace"));
  }

  // Register time derivatives
  TimeDomainProblemBuilder::RegisterGridFunctions();
};

void
DualFormulation::RegisterCoefficients()
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

WeakCurlEquationSystem::WeakCurlEquationSystem(const hephaestus::InputParameters & params)
  : TimeDependentEquationSystem(params),
    _h_curl_var_name(params.GetParam<std::string>("HCurlVarName")),
    _h_div_var_name(params.GetParam<std::string>("HDivVarName")),
    _alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
    _beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
    _dtalpha_coef_name(std::string("dt_") + _alpha_coef_name)
{
}

void
WeakCurlEquationSystem::Init(hephaestus::GridFunctions & gridfunctions,
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
WeakCurlEquationSystem::AddKernels()
{
  AddVariableNameIfMissing(_h_curl_var_name);
  AddVariableNameIfMissing(_h_div_var_name);
  std::string dh_curl_var_dt = GetTimeDerivativeName(_h_curl_var_name);

  // (αv_{n}, ∇×u')
  hephaestus::InputParameters weak_curl_params;
  weak_curl_params.SetParam("HCurlVarName", _h_curl_var_name);
  weak_curl_params.SetParam("HDivVarName", _h_div_var_name);
  weak_curl_params.SetParam("CoefficientName", _alpha_coef_name);
  AddKernel(_h_curl_var_name, std::make_unique<hephaestus::WeakCurlKernel>(weak_curl_params));

  // (αdt∇×u_{n+1}, ∇×u')
  hephaestus::InputParameters curl_curl_params;
  curl_curl_params.SetParam("CoefficientName", _dtalpha_coef_name);
  AddKernel(_h_curl_var_name, std::make_unique<hephaestus::CurlCurlKernel>(curl_curl_params));

  // (βu_{n+1}, u')
  hephaestus::InputParameters vector_fe_mass_params;
  vector_fe_mass_params.SetParam("CoefficientName", _beta_coef_name);
  AddKernel(_h_curl_var_name,
            std::make_unique<hephaestus::VectorFEMassKernel>(vector_fe_mass_params));
}

DualOperator::DualOperator(mfem::ParMesh & pmesh,
                           hephaestus::FESpaces & fespaces,
                           hephaestus::GridFunctions & gridfunctions,
                           hephaestus::BCMap & bc_map,
                           hephaestus::Coefficients & coefficients,
                           hephaestus::Sources & sources,
                           hephaestus::InputParameters & solver_options)
  : TimeDomainEquationSystemOperator(
        pmesh, fespaces, gridfunctions, bc_map, coefficients, sources, solver_options)
{
}

void
DualOperator::Init(mfem::Vector & X)
{
  TimeDomainEquationSystemOperator::Init(X);
  auto * eqs = dynamic_cast<hephaestus::WeakCurlEquationSystem *>(_equation_system);
  _h_curl_var_name = eqs->_h_curl_var_name;
  _h_div_var_name = eqs->_h_div_var_name;
  _u = _gridfunctions.Get(_h_curl_var_name);
  _dv = _gridfunctions.Get(GetTimeDerivativeName(_h_div_var_name));
  _h_curl_fe_space = _u->ParFESpace();
  _h_div_fe_space = _dv->ParFESpace();

  _curl = std::make_unique<mfem::ParDiscreteLinearOperator>(_h_curl_fe_space, _h_div_fe_space);
  _curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  _curl->Assemble();
}

void
DualOperator::ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt)
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

  _jacobian_solver =
      std::make_unique<hephaestus::DefaultHCurlPCGSolver>(_solver_options,
                                                          *_block_a.As<mfem::HypreParMatrix>(),
                                                          _equation_system->_test_pfespaces.at(0));

  _jacobian_solver->Mult(_true_rhs, _true_x);
  _equation_system->RecoverFEMSolution(_true_x, _gridfunctions);

  // Subtract off contribution from source
  _sources.SubtractSources(_u);

  // dv/dt_{n+1} = -∇×u
  _curl->Mult(*_u, *_dv);
  *_dv *= -1.0;
}

void
DualOperator::SetGridFunctions()
{
  TimeDomainEquationSystemOperator::SetGridFunctions();
  // Blocks for solution vector are smaller than the operator size
  // for DualOperator, as curl is stored separately.
  // Block operator only has the HCurl TrueVSize;
  _block_true_offsets.SetSize(_local_test_vars.size());
  _block_true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _local_test_vars.size() - 1; ++ind)
  {
    _block_true_offsets[ind + 1] = _local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  _block_true_offsets.PartialSum();

  _true_x.Update(_block_true_offsets);
  _true_rhs.Update(_block_true_offsets);
}

} // namespace hephaestus

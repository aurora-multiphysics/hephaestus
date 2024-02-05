#include "time_domain_problem_builder.hpp"

namespace hephaestus
{

TimeDomainProblem::~TimeDomainProblem()
{
  _td_equation_system.reset();
  _td_operator.reset();
}

std::vector<mfem::ParGridFunction *>
TimeDomainProblemBuilder::RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                                                  hephaestus::GridFunctions & gridfunctions)
{
  std::vector<mfem::ParGridFunction *> time_derivatives;

  for (auto & gridfunction_name : gridfunction_names)
  {
    gridfunctions.Register(GetTimeDerivativeName(gridfunction_name),
                           std::make_shared<mfem::ParGridFunction>(
                               gridfunctions.GetPtr(gridfunction_name, false)->ParFESpace()));

    time_derivatives.push_back(
        gridfunctions.GetPtr(GetTimeDerivativeName(gridfunction_name), false));
  }

  return time_derivatives;
}

void
TimeDomainProblemBuilder::RegisterGridFunctions()
{
  std::vector<std::string> gridfunction_names;
  for (auto const & [name, gf] : _problem->_gridfunctions)
  {
    gridfunction_names.push_back(name);
  }
  RegisterTimeDerivatives(gridfunction_names, _problem->_gridfunctions);
}

void
TimeDomainProblemBuilder::ConstructEquationSystem()
{
  hephaestus::InputParameters params;
  _problem->_td_equation_system = std::make_unique<hephaestus::TimeDependentEquationSystem>(params);
}

void
TimeDomainProblemBuilder::InitializeKernels()
{
  _problem->_td_equation_system->Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);

  _problem->_preprocessors.Init(_problem->_gridfunctions, _problem->_coefficients);
  _problem->_sources.Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
}

void
TimeDomainProblemBuilder::ConstructOperator()
{
  _problem->_td_operator =
      std::make_unique<hephaestus::TimeDomainEquationSystemOperator>(*(_problem->_pmesh),
                                                                     _problem->_fespaces,
                                                                     _problem->_gridfunctions,
                                                                     _problem->_bc_map,
                                                                     _problem->_coefficients,
                                                                     _problem->_sources,
                                                                     _problem->_solver_options);
  _problem->_td_operator->SetEquationSystem(_problem->_td_equation_system.get());
  _problem->_td_operator->SetGridFunctions();
}

void
TimeDomainProblemBuilder::ConstructState()
{
  // Vector of dofs.
  _problem->_f = std::make_unique<mfem::BlockVector>(_problem->_td_operator->_true_offsets);
  *(_problem->_f) = 0.0;                         // give initial value
  _problem->_td_operator->Init(*(_problem->_f)); // Set up initial conditions
  _problem->_td_operator->SetTime(0.0);
}

void
TimeDomainProblemBuilder::ConstructSolver()
{
  _problem->_ode_solver = std::make_unique<mfem::BackwardEulerSolver>();
  _problem->_ode_solver->Init(*(_problem->_td_operator));
}

} // namespace hephaestus

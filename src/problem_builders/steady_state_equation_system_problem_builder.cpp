#include "steady_state_equation_system_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateEquationSystemProblemBuilder::SetOperatorGridFunctions()
{
  _problem->GetOperator()->SetGridFunctions();
}

void
SteadyStateEquationSystemProblemBuilder::InitializeKernels()
{
  _problem->GetEquationSystem()->Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
  _problem->_preprocessors.Init(_problem->_gridfunctions, _problem->_coefficients);
  _problem->_sources.Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
}

void
SteadyStateEquationSystemProblemBuilder::ConstructOperator()
{
  _problem->SetOperator(std::make_unique<hephaestus::EquationSystemProblemOperator>(*_problem));

  auto equation_system = std::make_unique<hephaestus::EquationSystem>();

  _problem->GetOperator()->SetEquationSystem(std::move(equation_system));
}

void
SteadyStateEquationSystemProblemBuilder::ConstructState()
{
  _problem->_f =
      std::make_unique<mfem::BlockVector>(_problem->GetOperator()->_true_offsets); // Vector of dofs
  _problem->GetOperator()->Init(*(_problem->_f)); // Set up initial conditions
}

} // namespace hephaestus

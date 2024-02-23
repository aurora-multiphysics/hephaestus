#include "steady_state_equation_system_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateEquationSystemProblemBuilder::SetOperatorGridFunctions()
{
  _problem->GetOperator()->SetGridFunctions();
}

void
SteadyStateEquationSystemProblemBuilder::ConstructOperator()
{
  auto equation_system = std::make_unique<hephaestus::EquationSystem>();

  _problem->SetOperator(std::make_unique<hephaestus::EquationSystemProblemOperator>(
      *_problem, std::move(equation_system)));
}

void
SteadyStateEquationSystemProblemBuilder::ConstructState()
{
  _problem->_f =
      std::make_unique<mfem::BlockVector>(_problem->GetOperator()->_true_offsets); // Vector of dofs
  _problem->GetOperator()->Init(*(_problem->_f)); // Set up initial conditions
}

} // namespace hephaestus

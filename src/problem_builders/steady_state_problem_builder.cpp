#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::SetOperatorGridFunctions()
{
  _problem->GetOperator()->SetGridFunctions();
}

void
SteadyStateProblemBuilder::ConstructOperator()
{
  _problem->ConstructOperator();
}

void
SteadyStateProblemBuilder::ConstructState()
{
  _problem->_f =
      std::make_unique<mfem::BlockVector>(_problem->GetOperator()->_true_offsets); // Vector of dofs
  _problem->GetOperator()->Init(*(_problem->_f)); // Set up initial conditions
}
} // namespace hephaestus

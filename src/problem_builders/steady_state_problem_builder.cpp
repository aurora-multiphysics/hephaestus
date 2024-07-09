#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::ConstructOperator()
{
  GetProblem()->ConstructOperator();
}

} // namespace hephaestus

#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::ConstructOperator()
{
  auto problem_operator = std::make_unique<hephaestus::ProblemOperator>(*GetProblem());
  GetProblem()->SetOperator(std::move(problem_operator));
}

} // namespace hephaestus

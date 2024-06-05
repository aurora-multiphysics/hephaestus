#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::ConstructOperator()
{
  InputParameters params;
  params.Set("Problem", GetBaseProblem());

  auto problem_operator = std::make_unique<hephaestus::ProblemOperator>(params);
  GetProblem()->SetOperator(std::move(problem_operator));
}

} // namespace hephaestus

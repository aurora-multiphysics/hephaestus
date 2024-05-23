#include "steady_state_equation_system_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateEquationSystemProblemBuilder::ConstructOperator()
{
  auto equation_system = std::make_unique<hephaestus::EquationSystem>();
  auto problem_operator = std::make_unique<hephaestus::EquationSystemProblemOperator>(
      *GetProblem(), std::move(equation_system));

  GetProblem()->SetOperator(std::move(problem_operator));
}

} // namespace hephaestus

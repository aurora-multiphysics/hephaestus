#include "time_domain_equation_system_problem_builder.hpp"

namespace hephaestus
{
void
TimeDomainEquationSystemProblemBuilder::ConstructOperator()
{
  auto equation_system = std::make_unique<hephaestus::TimeDependentEquationSystem>();
  auto problem_operator = std::make_unique<hephaestus::TimeDomainEquationSystemProblemOperator>(
      *GetProblem(), std::move(equation_system));

  GetProblem()->SetOperator(std::move(problem_operator));
}

} // namespace hephaestus

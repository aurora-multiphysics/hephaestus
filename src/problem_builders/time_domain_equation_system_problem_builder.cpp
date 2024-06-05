#include "time_domain_equation_system_problem_builder.hpp"

namespace hephaestus
{
void
TimeDomainEquationSystemProblemBuilder::ConstructOperator()
{
  InputParameters params;
  params.Set("Problem", GetBaseProblem());

  auto equation_system = std::make_unique<hephaestus::TimeDependentEquationSystem>();
  auto problem_operator = std::make_unique<hephaestus::TimeDomainEquationSystemProblemOperator>(
      params, std::move(equation_system));

  GetProblem()->SetOperator(std::move(problem_operator));
}

} // namespace hephaestus

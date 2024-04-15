#include "equation_system_problem_operator.hpp"

namespace hephaestus
{
void
EquationSystemProblemOperator::SetGridFunctions()
{
  _trial_var_names = GetEquationSystem()->_trial_var_names;
  ProblemOperator::SetGridFunctions();
}
}
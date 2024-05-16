#include "equation_system_problem_operator.hpp"

namespace hephaestus
{
void
EquationSystemProblemOperator::SetTrialVariableNames()
{
  _trial_var_names = GetEquationSystem()->_trial_var_names;
}

void
EquationSystemProblemOperator::Init()
{
  GetEquationSystem()->Init(_problem._gridfunctions,
                            _problem._fespaces,
                            _problem._bc_map,
                            _problem._coefficients,
                            _problem._sources);

  ProblemOperator::Init();
}

} // namespace hephaestus
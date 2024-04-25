#include "equation_system_problem_operator.hpp"

namespace hephaestus
{
void
EquationSystemProblemOperator::SetTrialVariables()
{
  _trial_var_names = GetEquationSystem()->_trial_var_names;
  ProblemOperator::SetTrialVariables();
}

void
EquationSystemProblemOperator::Init(mfem::BlockVector & X)
{
  ProblemOperator::Init(X);

  GetEquationSystem()->BuildEquationSystem(_problem._bc_map, _problem._sources);
}

void
EquationSystemProblemOperator::Init()
{
  GetEquationSystem()->Init(
      _problem._gridfunctions, _problem._fespaces, _problem._bc_map, _problem._coefficients);

  ProblemOperator::Init();
}

} // namespace hephaestus
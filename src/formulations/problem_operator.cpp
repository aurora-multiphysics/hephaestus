#include "problem_operator.hpp"

namespace hephaestus
{

void
ProblemOperator::SetGridFunctions()
{
  _trial_variables = _problem._gridfunctions.Get(_trial_var_names);

  // Set operator size and block structure
  _block_true_offsets.SetSize(_trial_variables.size() + 1);
  _block_true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _block_true_offsets[ind + 1] = _trial_variables.at(ind)->ParFESpace()->TrueVSize();
  }
  _block_true_offsets.PartialSum();

  _true_offsets.SetSize(_trial_variables.size() + 1);
  _true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _true_offsets[ind + 1] = _trial_variables.at(ind)->ParFESpace()->GetVSize();
  }
  _true_offsets.PartialSum();

  height = _true_offsets[_trial_variables.size()];
  width = _true_offsets[_trial_variables.size()];
  _true_x.Update(_block_true_offsets);
  _true_rhs.Update(_block_true_offsets);
};

void
ProblemOperator::Init(mfem::Vector & X)
{
  // Define material property coefficients
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _trial_variables.at(ind)->MakeRef(
        _trial_variables.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
  }
}

EquationSystemProblemOperator::EquationSystemProblemOperator(hephaestus::Problem & problem)
  : ProblemOperator(problem)
{
}

} // namespace hephaestus
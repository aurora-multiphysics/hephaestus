#include "problem_operator.hpp"

namespace hephaestus
{

void
ProblemOperator::SetGridFunctions()
{
  ProblemOperatorInterface::SetGridFunctions();
  width = height = _true_offsets[_trial_variables.size()];
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
} // namespace hephaestus
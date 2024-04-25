#include "problem_operator_interface.hpp"

namespace hephaestus
{

void
ProblemOperatorInterface::SetTrialVariables()
{
  _trial_variables = _problem._gridfunctions.Get(_trial_var_names);
}

void
ProblemOperatorInterface::UpdateOffsets()
{
  return UpdateOffsetsWithSize(_trial_variables.size());
}

void
ProblemOperatorInterface::UpdateOffsetsWithSize(const size_t soln_vector_size)
{
  if (soln_vector_size > _trial_variables.size())
  {
    MFEM_ABORT("Solution vector size (" << soln_vector_size
                                        << ") cannot exceed the number of trial variables ("
                                        << _trial_variables.size() << ").");
  }

  _block_true_offsets.SetSize(soln_vector_size + 1);
  _true_offsets.SetSize(soln_vector_size + 1);

  _block_true_offsets[0] = _true_offsets[0] = 0;
  for (size_t i = 0; i < soln_vector_size; ++i)
  {
    mfem::ParFiniteElementSpace * fespace = _trial_variables.at(i)->ParFESpace();

    _block_true_offsets[i + 1] = fespace->TrueVSize();
    _true_offsets[i + 1] = fespace->GetVSize();
  }

  // Partial sum over values to calculate offsets.
  _block_true_offsets.PartialSum();
  _true_offsets.PartialSum();
}

void
ProblemOperatorInterface::UpdateBlockVectors()
{
  _true_x.Update(_block_true_offsets);
  _true_rhs.Update(_block_true_offsets);
}

void
ProblemOperatorInterface::SetGridFunctions()
{
  SetTrialVariables();

  // Recalculate the offsets from gridfunction trial variables.
  UpdateOffsets();

  // Update the block vectors with new offsets.
  UpdateBlockVectors();
}

void
ProblemOperatorInterface::Init(mfem::BlockVector & X)
{
  X.Update(_true_offsets);

  for (size_t i = 0; i < _trial_variables.size(); ++i)
  {
    mfem::ParGridFunction * trial_var = _trial_variables.at(i);

    trial_var->MakeRef(trial_var->ParFESpace(), X, _true_offsets[i]);
    *trial_var = 0.0;
  }
}

void
ProblemOperatorInterface::Init()
{
  SetGridFunctions();
  Init(*_problem._f);
}

void
ProblemOperatorInterface::Update()
{
  logger.debug("Update called for ProblemOperatorInterface.");
  Init();
}

}
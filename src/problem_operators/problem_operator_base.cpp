#include "problem_operator_base.hpp"

namespace hephaestus
{

void
ProblemOperatorBase::SetTrialVariables()
{
  SetTrialVariableNames();

  _trial_variables = _problem._gridfunctions.Get(_trial_var_names);
}

void
ProblemOperatorBase::UpdateOffsets()
{
  UpdateOffsetsWithSize(_trial_variables.size());
}

void
ProblemOperatorBase::UpdateOffsetsWithSize(const size_t block_vector_size)
{
  if (block_vector_size > _trial_variables.size())
  {
    MFEM_ABORT("Solution vector size (" << block_vector_size
                                        << ") cannot exceed the number of trial variables ("
                                        << _trial_variables.size() << ").");
  }

  _block_true_offsets.SetSize(block_vector_size + 1);
  _true_offsets.SetSize(_trial_variables.size() + 1);

  _block_true_offsets[0] = _true_offsets[0] = 0;
  for (size_t i = 0; i < _trial_variables.size(); ++i)
  {
    mfem::ParFiniteElementSpace * fespace = _trial_variables.at(i)->ParFESpace();

    if (i < block_vector_size)
      _block_true_offsets[i + 1] = fespace->GetTrueVSize();

    _true_offsets[i + 1] = fespace->GetVSize();
  }

  // Partial sum over values to calculate offsets.
  _block_true_offsets.PartialSum();
  _true_offsets.PartialSum();

  // Update block vectors.
  _true_x.Update(_block_true_offsets);
  _true_rhs.Update(_block_true_offsets);

  // Set width and height member variables.
  Width() = Height() = _true_offsets.Last();
}

void
ProblemOperatorBase::UpdateBlockVector(mfem::BlockVector & X)
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
ProblemOperatorBase::Init()
{
  SetTrialVariables();

  UpdateOffsets();
  UpdateBlockVector(*_problem._f);
}

void
ProblemOperatorBase::Update()
{
  logger.debug("Update called for ProblemOperatorBase.");

  // Recalculate the offsets from gridfunction trial variables.
  UpdateOffsets();

  UpdateBlockVector(*_problem._f);
}

}
#include "problem_operator_base.hpp"

namespace hephaestus
{

ProblemOperatorBase::ProblemOperatorBase(hephaestus::Problem & problem) : _problem{problem}
{
  _block_vector = std::make_unique<mfem::BlockVector>();
}

void
ProblemOperatorBase::SetTrialVariables()
{
  SetTrialVariableNames();

  _trial_variables = _problem._gridfunctions.Get(_trial_var_names);
}

int
ProblemOperatorBase::GetSolutionVectorSize() const
{
  return static_cast<int>(GetTrialVariablesSize());
}

int
ProblemOperatorBase::GetTrialVariablesSize() const
{
  return static_cast<int>(_trial_variables.size());
}

void
ProblemOperatorBase::UpdateOffsets()
{
  if (GetSolutionVectorSize() > GetTrialVariablesSize())
  {
    MFEM_ABORT("Solution vector size (" << GetSolutionVectorSize()
                                        << ") cannot exceed the number of trial variables ("
                                        << GetTrialVariablesSize() << ").");
  }

  _block_true_offsets.SetSize(GetSolutionVectorSize() + 1);
  _true_offsets.SetSize(GetTrialVariablesSize() + 1);

  _block_true_offsets[0] = _true_offsets[0] = 0;
  for (int i = 0; i < GetTrialVariablesSize(); ++i)
  {
    mfem::ParFiniteElementSpace * fespace = _trial_variables.at(i)->ParFESpace();

    if (i < GetSolutionVectorSize())
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

  for (int i = 0; i < GetTrialVariablesSize(); ++i)
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
  UpdateBlockVector(*_block_vector);
}

void
ProblemOperatorBase::Update()
{
  logger.debug("Update called for ProblemOperatorBase.");

  // Recalculate the offsets from gridfunction trial variables.
  UpdateOffsets();

  UpdateBlockVector(*_block_vector);
}

}
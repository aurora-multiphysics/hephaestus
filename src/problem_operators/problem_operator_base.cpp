#include "problem_operator_base.hpp"

namespace hephaestus
{

ProblemOperatorBase::ProblemOperatorBase(hephaestus::Problem & problem) : _problem{problem} {}

void
ProblemOperatorBase::ConstructJacobianSolver()
{
  auto precond = std::make_unique<mfem::HypreBoomerAMG>();
  precond->SetPrintLevel(GetGlobalPrintLevel());

  auto solver = std::make_unique<mfem::HypreGMRES>(_problem._comm);

  solver->SetTol(1e-16);
  solver->SetAbsTol(1e-16);
  solver->SetMaxIter(1000);
  solver->SetKDim(10);
  solver->SetPrintLevel(GetGlobalPrintLevel());
  solver->SetPreconditioner(*precond);

  _jacobian_preconditioner = std::move(precond);
  _jacobian_solver = std::move(solver);
}

void
ProblemOperatorBase::ConstructNonlinearSolver()
{
  _nonlinear_solver = std::make_unique<mfem::NewtonSolver>(_problem._comm);

  // Defaults to one iteration, without further nonlinear iterations
  _nonlinear_solver->SetRelTol(0.0);
  _nonlinear_solver->SetAbsTol(0.0);
  _nonlinear_solver->SetMaxIter(1);
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
  _block_vector = std::make_unique<mfem::BlockVector>();
  ConstructJacobianSolver();
  ConstructNonlinearSolver();

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
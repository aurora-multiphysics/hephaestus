#include "problem_operator.hpp"

namespace hephaestus
{

ProblemOperator::ProblemOperator(hephaestus::Problem & problem) : _problem{problem} {}

ProblemOperator::SolverOptions
ProblemOperator::DefaultSolverOptions() const
{
  return {._tolerance = 1e-16,
          ._abs_tolerance = 1e-16,
          ._max_iteration = 1000,
          ._print_level = GetGlobalPrintLevel(),
          ._k_dim = 10};
}

void
ProblemOperator::Solve()
{
  MFEM_ABORT("Not implemented.");
}

void
ProblemOperator::Mult(const mfem::Vector & x, mfem::Vector & y) const
{
  MFEM_ABORT("Not implemented.");
}

void
ProblemOperator::ConstructJacobianSolver()
{
  auto precond = std::make_unique<mfem::HypreBoomerAMG>();
  precond->SetPrintLevel(GetGlobalPrintLevel());

  auto solver = std::make_unique<mfem::HypreGMRES>(_problem._comm);
  solver->SetPreconditioner(*precond);

  _jacobian_preconditioner = std::move(precond);
  _jacobian_solver = std::move(solver);
}

void
ProblemOperator::ConstructJacobianSolverAndApplyOptions()
{
  ConstructJacobianSolver();
  ApplySolverOptions();
}

void
ProblemOperator::SetSolverOptions(SolverOptions options)
{
  // Store the options for future.
  _solver_options = std::move(options);
  ApplySolverOptions();
}

void
ProblemOperator::ApplySolverOptions()
{
  auto & solver = static_cast<mfem::HypreGMRES &>(*_jacobian_solver);

  solver.SetTol(GetSolverOptions()._tolerance);
  solver.SetAbsTol(GetSolverOptions()._abs_tolerance);
  solver.SetMaxIter(GetSolverOptions()._max_iteration);
  solver.SetKDim(GetSolverOptions()._k_dim);
  solver.SetPrintLevel(GetSolverOptions()._print_level);
}

void
ProblemOperator::ConstructNonlinearSolver()
{
  _nonlinear_solver = std::make_unique<mfem::NewtonSolver>(_problem._comm);

  // Defaults to one iteration, without further nonlinear iterations
  _nonlinear_solver->SetRelTol(0.0);
  _nonlinear_solver->SetAbsTol(0.0);
  _nonlinear_solver->SetMaxIter(1);
}

void
ProblemOperator::SetTrialVariables()
{
  SetTrialVariableNames();

  _trial_variables = _problem._gridfunctions.Get(_trial_var_names);
}

int
ProblemOperator::GetSolutionVectorSize() const
{
  return static_cast<int>(GetTrialVariablesSize());
}

int
ProblemOperator::GetTrialVariablesSize() const
{
  return static_cast<int>(_trial_variables.size());
}

void
ProblemOperator::UpdateOffsets()
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
ProblemOperator::UpdateBlockVector(mfem::BlockVector & X)
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
ProblemOperator::Init()
{
  _block_vector = std::make_unique<mfem::BlockVector>();
  _solver_options = DefaultSolverOptions();

  ConstructJacobianSolverAndApplyOptions();
  ConstructNonlinearSolver();

  SetTrialVariables();

  UpdateOffsets();
  UpdateBlockVector(*_block_vector);
}

void
ProblemOperator::Update()
{
  logger.debug("Update called for ProblemOperator.");

  // Recalculate the offsets from gridfunction trial variables.
  UpdateOffsets();

  UpdateBlockVector(*_block_vector);

  // Update the Jacobian solver.
  ConstructJacobianSolverAndApplyOptions();
}

}
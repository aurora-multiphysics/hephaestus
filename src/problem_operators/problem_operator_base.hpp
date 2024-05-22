#pragma once
#include "mfem.hpp"
#include "problem_builder_base.hpp"
#include <vector>
#include <string>

namespace hephaestus
{
class ProblemOperatorBase
{
public:
  ProblemOperatorBase() = delete;

  virtual ~ProblemOperatorBase() = default;

  /// Initialize the problem operator.
  virtual void Init();

  /// Update the problem operator after a mesh change.
  virtual void Update();

  /// Set the Jacobian preconditioner.
  void SetJacobianPreconditioner(std::unique_ptr<mfem::Solver> preconditioner)
  {
    _jacobian_preconditioner = std::move(preconditioner);
  }

  /// Set the Jacobian solver.
  void SetJacobianSolver(std::unique_ptr<mfem::Solver> solver)
  {
    _jacobian_solver = std::move(solver);
  }

  /// Set the nonlinear solver.
  void SetNonlinearSolver(std::unique_ptr<mfem::NewtonSolver> nl_solver)
  {
    _nonlinear_solver = std::move(nl_solver);
  }

  /// Accessor for Jacobian preconditioner.
  template <class TSolver>
  TSolver * JacobianPreconditioner() const
  {
    return static_cast<TSolver *>(_jacobian_preconditioner.get());
  }

protected:
  /// Use of protected constructor to only allow construction by derived classes.
  /// All problem operator classes are built on-top of this class and it should not
  /// be possible to use directly.
  explicit ProblemOperatorBase(hephaestus::Problem & problem);

  /// Set trial variables names. Override in derived classes.
  virtual void SetTrialVariableNames() {}

  /// Set gridfunction of trial variables from the trial variable names.
  virtual void SetTrialVariables();

  /// Returns a reference to the operator's width.
  virtual int & Width() = 0;

  /// Returns a reference to the operator's height.
  virtual int & Height() = 0;

  /// Returns the solution vector size. By default this will be the same as the
  /// number of trial variables.
  virtual int GetSolutionVectorSize() const;

  /// Returns the number of trial variables.
  int GetTrialVariablesSize() const;

  // Reference to the current problem.
  hephaestus::Problem & _problem;

  /// Vector of names of state gridfunctions used in formulation, ordered by
  /// appearance in block vector during solve.
  std::vector<std::string> _trial_var_names;

  /// Vector of trial variables set by SetTrialVariables.
  std::vector<mfem::ParGridFunction *> _trial_variables;

  /// Arrays of offsets.
  mfem::Array<int> _true_offsets, _block_true_offsets;

  /// Block vectors.
  mfem::BlockVector _true_x, _true_rhs;

  /// Store the Jacobian preconditioner.
  std::unique_ptr<mfem::Solver> _jacobian_preconditioner{nullptr};

  /// Store the Jacobian solver.
  std::unique_ptr<mfem::Solver> _jacobian_solver{nullptr};

  /// Store the non-linear solver.
  std::unique_ptr<mfem::NewtonSolver> _nonlinear_solver{nullptr};

  /// Store the block vector.
  std::unique_ptr<mfem::BlockVector> _block_vector{nullptr};

private:
  /// Update the block vectors and offsets after a mesh change.
  void UpdateOffsets();

  /// Update a block vector. Should be called after the offsets have been updated.
  void UpdateBlockVector(mfem::BlockVector & X);
};
}
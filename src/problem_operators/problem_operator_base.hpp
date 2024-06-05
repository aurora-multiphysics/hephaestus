#pragma once
#include "mfem.hpp"
#include "problem_builder_base.hpp"
#include <vector>
#include <string>

namespace hephaestus
{
/// Steady-State Problem Operator with no equation system.
class ProblemOperator : virtual public mfem::Operator
{
public:
  // ProblemOperator() = delete;
  explicit ProblemOperator(hephaestus::Problem & problem);

  ~ProblemOperator() override = default;

  /// Initialize the problem operator.
  virtual void Init();

  /// Update the problem operator after a mesh change.
  virtual void Update();

  /// A structure for setting solver options.
  struct SolverOptions
  {
    double _tolerance;
    double _abs_tolerance;

    unsigned int _max_iteration;

    int _print_level;
    int _k_dim;
  };

  /// Sets the solver's options. Then calls ApplySolverOptions.
  void SetSolverOptions(SolverOptions options);

  virtual void Solve() { Solve(*_block_vector); }

  virtual void Solve(mfem::Vector & X) {}

  void Mult(const mfem::Vector & x, mfem::Vector & y) const override {}

protected:
  /// Use of protected constructor to only allow construction by derived classes.
  /// All problem operator classes are built on-top of this class and it should not
  /// be possible to use directly.
  // explicit ProblemOperator(hephaestus::Problem & problem);

  /// Override in derived classes to set default solver options. These will be
  /// applied following solver construction.
  virtual SolverOptions DefaultSolverOptions() const;

  /// Applies the current solver options. Override in derived classes.
  virtual void ApplySolverOptions();

  /// Override in derived classes to construct the Jacobian solver. Called in the
  /// Init method.
  virtual void ConstructJacobianSolver();

  /// Override in derived classes to construct the non-linear solver. Called in
  /// the Init method.
  virtual void ConstructNonlinearSolver();

  /// Set trial variables names. Override in derived classes.
  virtual void SetTrialVariableNames() {}

  /// Set gridfunction of trial variables from the trial variable names.
  virtual void SetTrialVariables();

  /// Returns a reference to the operator's width.
  virtual int & Width() { return mfem::Operator::width; }

  /// Returns a reference to the operator's height.
  virtual int & Height() { return mfem::Operator::height; }

  /// Returns the solution vector size. By default this will be the same as the
  /// number of trial variables.
  virtual int GetSolutionVectorSize() const;

  /// Returns the number of trial variables.
  int GetTrialVariablesSize() const;

  /// Solver options accessor.
  const SolverOptions & GetSolverOptions() const { return _solver_options; }

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
  /// Calls ConstructJacobianSolver followed by ApplySolverOptions. This ensures
  /// that ApplySolverOptions is always called!
  void ConstructJacobianSolverAndApplyOptions();

  /// Update the block vectors and offsets after a mesh change.
  void UpdateOffsets();

  /// Update a block vector. Should be called after the offsets have been updated.
  void UpdateBlockVector(mfem::BlockVector & X);

  /// The current solver options.
  SolverOptions _solver_options;
};
}
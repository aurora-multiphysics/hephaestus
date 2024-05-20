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

protected:
  /// Use of protected constructor to only allow construction by derived classes.
  /// All problem operator classes are built on-top of this class and it should not
  /// be possible to use directly.
  explicit ProblemOperatorBase(hephaestus::Problem & problem) : _problem(problem) {}

  /// Set trial variables names. Override in derived classes.
  virtual void SetTrialVariableNames() {}

  /// Set gridfunction of trial variables from the trial variable names.
  virtual void SetTrialVariables();

  /// Update the block vectors and offsets after a mesh change.
  virtual void UpdateOffsets();

  /// Update a block vector. Should be called after the offsets have been updated.
  virtual void UpdateBlockVector(mfem::BlockVector & X);

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
};
}
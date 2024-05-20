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
  ProblemOperatorBase(hephaestus::Problem & problem) : _problem(problem) {}

  virtual ~ProblemOperatorBase() = default;

  /// Initialize the problem operator.
  virtual void Init();

  /// Update the problem operator after a mesh change.
  virtual void Update();

  /// Arrays of offsets.
  mfem::Array<int> _true_offsets, _block_true_offsets;

  /// Block vectors.
  mfem::BlockVector _true_x, _true_rhs;

protected:
  /// Set trial variables names. Override in derived classes.
  virtual void SetTrialVariableNames() {}

  /// Set gridfunction of trial variables from the trial variable names.
  virtual void SetTrialVariables();

  /// Update the block vectors and offsets after a mesh change.
  virtual void UpdateOffsets();

  /// Update a block vector. Should be called after the offsets have been updated.
  virtual void UpdateBlockVector(mfem::BlockVector & X);

  /// Called by UpdateOffsets.
  virtual void UpdateOffsetsWithSize(size_t soln_vector_size);

  /// Returns a reference to the operator's width.
  virtual int & Width() = 0;

  /// Returns a reference to the operator's height.
  virtual int & Height() = 0;

  // Reference to the current problem.
  hephaestus::Problem & _problem;

  /// Vector of names of state gridfunctions used in formulation, ordered by
  /// appearance in block vector during solve.
  std::vector<std::string> _trial_var_names;

  /// Vector of trial variables set by SetTrialVariables.
  std::vector<mfem::ParGridFunction *> _trial_variables;
};
}
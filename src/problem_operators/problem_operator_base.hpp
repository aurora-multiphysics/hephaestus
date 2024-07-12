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

  virtual void Init();

  virtual void Update();

  mfem::Array<int> _true_offsets, _block_true_offsets;

  mfem::BlockVector _true_x, _true_rhs;

protected:
  virtual void SetTrialVariableNames() {}
  virtual void SetTrialVariables();
  virtual void UpdateOffsets();

  virtual void UpdateBlockVector(mfem::BlockVector & X);

  virtual void UpdateOffsetsWithSize(size_t soln_vector_size);

  /// Returns a reference to the operator's width.
  virtual int & Width() = 0;

  /// Returns a reference to the operator's height.
  virtual int & Height() = 0;

  // Reference to the current problem.
  hephaestus::Problem & _problem;

  // Vector of names of state gridfunctions used in formulation, ordered by appearance in block
  // vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;
};
}
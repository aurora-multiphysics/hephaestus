#pragma once
#include "problem_operator_interface.hpp"

namespace hephaestus
{
class ProblemOperatorBase : public ProblemOperatorInterface
{
public:
  ProblemOperatorBase(hephaestus::Problem & problem) : _problem(problem) {}
  ~ProblemOperatorBase() override = default;

  void UpdateOperatorWidthAndHeight() override {}

  void Init() override;

  void Update() override;

  mfem::Array<int> _true_offsets, _block_true_offsets;

  mfem::BlockVector _true_x, _true_rhs;

protected:
  void SetTrialVariableNames() override {}
  void SetTrialVariables() override;
  void UpdateOffsets() override;

  void UpdateBlockVector(mfem::BlockVector & X) override;

  void UpdateOffsetsWithSize(size_t soln_vector_size) override;

  // Reference to the current problem.
  hephaestus::Problem & _problem;

  // Vector of names of state gridfunctions used in formulation, ordered by appearance in block
  // vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;
};
}
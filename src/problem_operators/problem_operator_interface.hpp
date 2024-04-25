#pragma once
#include "problem_builder_base.hpp"

namespace hephaestus
{
/// Interface inherited by ProblemOperator and TimeDomainProblemOperator. Removes duplicated code in both classes.
class ProblemOperatorInterface
{
public:
  ProblemOperatorInterface(hephaestus::Problem & problem) : _problem(problem) {}
  virtual ~ProblemOperatorInterface() = default;

  virtual void Init(mfem::BlockVector & X);

  virtual void UpdateOperatorWidthAndHeight() {}

  virtual void Init();

  virtual void Update();

  mfem::Array<int> _true_offsets, _block_true_offsets;

  mfem::BlockVector _true_x, _true_rhs;
  mfem::OperatorHandle _equation_system_operator;

protected:
  virtual void SetTrialVariableNames() {}
  virtual void SetTrialVariables();
  virtual void UpdateOffsets();

  void UpdateOffsetsWithSize(size_t soln_vector_size);

  // Reference to the current problem.
  hephaestus::Problem & _problem;

  // Vector of names of state gridfunctions used in formulation, ordered by appearance in block
  // vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;
};
}
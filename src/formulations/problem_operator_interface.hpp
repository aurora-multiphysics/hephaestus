#pragma once
#include "problem_builder_base.hpp"

namespace hephaestus
{
class ProblemOperatorInterface
{
public:
  ProblemOperatorInterface(hephaestus::Problem & problem) : _problem(problem) {}
  virtual ~ProblemOperatorInterface() = default;

  [[nodiscard]] virtual mfem::Operator * GetEquationSystem() const = 0;

  virtual void SetGridFunctions() = 0;
  virtual void Init(mfem::Vector & X) = 0;

  mfem::Array<int> _true_offsets, _block_true_offsets;

  mfem::BlockVector _true_x, _true_rhs;
  mfem::OperatorHandle _equation_system_operator;

  // Vector of names of state gridfunctions used in formulation, ordered by appearance in block
  // vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;

protected:
  hephaestus::Problem & _problem;
};
}
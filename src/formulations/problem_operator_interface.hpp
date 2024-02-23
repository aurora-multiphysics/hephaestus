#pragma once
#include "problem_builder_base.hpp"

namespace hephaestus
{
/// Interface inherited by ProblemOperator and TimeDependentProblemOperator. Removes duplicated code in both classes.
class ProblemOperatorInterface
{
public:
  ProblemOperatorInterface(hephaestus::Problem & problem) : _problem(problem) {}
  virtual ~ProblemOperatorInterface() = default;

  virtual void SetGridFunctions() = 0;
  virtual void Init(mfem::Vector & X) = 0;

  mfem::Array<int> _true_offsets, _block_true_offsets;

  mfem::BlockVector _true_x, _true_rhs;
  mfem::OperatorHandle _equation_system_operator;

protected:
  hephaestus::Problem & _problem;

  // Vector of names of state gridfunctions used in formulation, ordered by appearance in block
  // vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;
};

/// Interface inherited by ProblemOperators with an equation system.
class EquationSystemProblemOperatorInterface : public ProblemOperatorInterface
{
public:
  EquationSystemProblemOperatorInterface(hephaestus::Problem & problem)
    : ProblemOperatorInterface(problem)
  {
  }

  ~EquationSystemProblemOperatorInterface() override = default;

  [[nodiscard]] virtual mfem::Operator * GetEquationSystem() const = 0;
};
}
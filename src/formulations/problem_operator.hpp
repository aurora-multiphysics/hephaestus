#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"

namespace hephaestus
{
class ProblemOperator : public mfem::Operator
{
public:
  ProblemOperator(hephaestus::Problem & problem) : _problem(problem){};

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector & X);
  virtual void Solve(mfem::Vector & X){};
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override{};

  mfem::Array<int> _true_offsets, _block_true_offsets;

  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;

  mfem::BlockVector _true_x, _true_rhs;
  mfem::OperatorHandle _equation_system_operator;

protected:
  hephaestus::Problem & _problem;
};

} // namespace hephaestus
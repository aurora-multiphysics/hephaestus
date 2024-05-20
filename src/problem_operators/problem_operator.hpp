#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"
#include "problem_operator_base.hpp"

namespace hephaestus
{
/// Steady-state problem operator with no equation system.
class ProblemOperator : public mfem::Operator, public ProblemOperatorBase
{
public:
  ProblemOperator(hephaestus::Problem & problem) : ProblemOperatorBase(problem) {}
  ~ProblemOperator() override = default;

  virtual void Solve(mfem::Vector & X) {}
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override {}

protected:
  int & Width() final { return mfem::Operator::width; }
  int & Height() final { return mfem::Operator::height; }
};

} // namespace hephaestus
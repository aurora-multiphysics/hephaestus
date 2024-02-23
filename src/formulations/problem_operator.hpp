#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"
#include "problem_operator_interface.hpp"

namespace hephaestus
{
/// Steady-state problem operator with no equation system.
class ProblemOperator : public mfem::Operator, public ProblemOperatorInterface
{
public:
  ProblemOperator(hephaestus::Problem & problem) : ProblemOperatorInterface(problem) {}
  ~ProblemOperator() override = default;

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;

  virtual void Solve(mfem::Vector & X) {}
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override {}
};

/// Steady-state problem operator with an equation system.
class EquationSystemProblemOperator : public ProblemOperator,
                                      public EquationSystemProblemOperatorInterface
{
public:
  EquationSystemProblemOperator(hephaestus::Problem &) = delete;
  EquationSystemProblemOperator(hephaestus::Problem & problem,
                                std::unique_ptr<EquationSystem> equation_system)
    : ProblemOperator(problem), _equation_system{std::move(equation_system)}
  {
  }

  ~EquationSystemProblemOperator() override = default;

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added to ProblemOperator.");
    }

    return _equation_system.get();
  }

private:
  std::unique_ptr<EquationSystem> _equation_system{nullptr};
};

} // namespace hephaestus
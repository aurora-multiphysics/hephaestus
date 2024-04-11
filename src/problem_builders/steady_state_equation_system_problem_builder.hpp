#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

/// Steady-state problems with an equation system.
class SteadyStateEquationSystemProblem : public SteadyStateProblem,
                                         public EquationSystemProblemInterface
{
public:
  SteadyStateEquationSystemProblem() = default;
  ~SteadyStateEquationSystemProblem() override = default;

  [[nodiscard]] EquationSystemProblemOperator * GetOperator() const override
  {
    return static_cast<EquationSystemProblemOperator *>(SteadyStateProblem::GetOperator());
  }

  void SetOperator(std::unique_ptr<EquationSystemProblemOperator> problem_operator)
  {
    SteadyStateProblem::SetOperator(std::move(problem_operator));
  }

  void ConstructOperator() override
  {
    auto equation_system = std::make_unique<EquationSystem>();
    auto problem_operator = std::make_unique<EquationSystemProblemOperator>(*this, std::move(equation_system));

    SetOperator(std::move(problem_operator));
  }

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }
};

/// Problem-builder for SteadyStateEquationSystemProblem.
class SteadyStateEquationSystemProblemBuilder : public SteadyStateProblemBuilder,
                                                public EquationSystemProblemBuilderInterface
{
public:
  /// NB: set "_problem" member variable in parent class.
  SteadyStateEquationSystemProblemBuilder() : SteadyStateProblemBuilder(std::make_unique<SteadyStateEquationSystemProblem>())
  {
  }

  ~SteadyStateEquationSystemProblemBuilder() override = default;

  /// NB: use of final! This calls ProblemBuilder::InitializeKernels and also ensures that the
  /// equation system is initialized.
  void InitializeKernels() final;

protected:
  [[nodiscard]] SteadyStateEquationSystemProblem * GetProblem() const override
  {
    return static_cast<SteadyStateEquationSystemProblem *>(SteadyStateProblemBuilder::GetProblem());
  }

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    return GetProblem()->GetEquationSystem();
  }
};

} // namespace hephaestus

#pragma once
#include "equation_system_problem_operator.hpp"
#include "problem_builder_base.hpp"
#include "steady_state_problem_builder.hpp"
#include "equation_system_interface.hpp"

namespace hephaestus
{

/// Steady-state problems with an equation system.
class SteadyStateEquationSystemProblem : public SteadyStateProblem, public EquationSystemInterface
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

  [[nodiscard]] hephaestus::EquationSystem * GetEquationSystem() const override
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
  SteadyStateEquationSystemProblemBuilder()
    : SteadyStateProblemBuilder(new SteadyStateEquationSystemProblem)
  {
  }

  ~SteadyStateEquationSystemProblemBuilder() override = default;

  auto ReturnProblem() { return ProblemBuilder::ReturnProblem<SteadyStateEquationSystemProblem>(); }

  void ConstructOperator() override;

protected:
  [[nodiscard]] hephaestus::SteadyStateEquationSystemProblem * GetProblem() const override
  {
    return ProblemBuilder::GetProblem<hephaestus::SteadyStateEquationSystemProblem>();
  }

  [[nodiscard]] hephaestus::EquationSystem * GetEquationSystem() const override
  {
    return GetProblem()->GetEquationSystem();
  }
};

} // namespace hephaestus

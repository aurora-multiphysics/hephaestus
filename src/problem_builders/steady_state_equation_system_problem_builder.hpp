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
    if (!_problem_operator)
    {
      MFEM_ABORT("No operator has been added.");
    }

    return _problem_operator.get();
  }

  void SetOperator(std::unique_ptr<EquationSystemProblemOperator> problem_operator)
  {
    _problem_operator.reset();
    _problem_operator = std::move(problem_operator);
  }

  void ConstructOperator() override
  {
    auto equation_system = std::make_unique<EquationSystem>();

    _problem_operator.reset();
    _problem_operator =
        std::make_unique<EquationSystemProblemOperator>(*this, std::move(equation_system));
  }

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }

private:
  std::unique_ptr<EquationSystemProblemOperator> _problem_operator{nullptr};
};

/// Problem-builder for SteadyStateEquationSystemProblem.
class SteadyStateEquationSystemProblemBuilder : public SteadyStateProblemBuilder,
                                                public EquationSystemProblemBuilderInterface
{
public:
  // NB: pass nullptr to parent class to prevent it constructing a private SteadyStateProblem
  // instance.
  SteadyStateEquationSystemProblemBuilder()
    : SteadyStateProblemBuilder(nullptr),
      _problem{std::make_unique<SteadyStateEquationSystemProblem>()}
  {
  }

  ~SteadyStateEquationSystemProblemBuilder() override = default;

  /// NB: use of final! This calls ProblemBuilder::InitializeKernels and also ensures that the
  /// equation system is initialized.
  void InitializeKernels() final;

  std::unique_ptr<SteadyStateEquationSystemProblem> ReturnProblem() { return std::move(_problem); }

protected:
  [[nodiscard]] SteadyStateEquationSystemProblem * GetProblem() const override
  {
    return _problem.get();
  }

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    return GetProblem()->GetEquationSystem();
  }

private:
  std::unique_ptr<SteadyStateEquationSystemProblem> _problem{nullptr};
};

} // namespace hephaestus

#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_builder.hpp"
#include "time_domain_problem_operator.hpp"

namespace hephaestus
{
class TimeDomainEquationSystemProblem : public TimeDomainProblem,
                                        public EquationSystemProblemInterface
{
public:
  [[nodiscard]] TimeDomainEquationSystemProblemOperator * GetOperator() const override
  {
    if (!_problem_operator)
    {
      MFEM_ABORT("No operator has been added.");
    }

    return _problem_operator.get();
  }

  void SetOperator(std::unique_ptr<TimeDomainEquationSystemProblemOperator> problem_operator)
  {
    _problem_operator.reset();
    _problem_operator = std::move(problem_operator);
  }

  void ConstructOperator() override
  {
    auto equation_system = std::make_unique<TimeDependentEquationSystem>();

    _problem_operator = std::make_unique<TimeDomainEquationSystemProblemOperator>(
        *this, std::move(equation_system));
  }

  [[nodiscard]] TimeDependentEquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }

private:
  std::unique_ptr<TimeDomainEquationSystemProblemOperator> _problem_operator{nullptr};
};

// Builder class of a time-domain EM formulation.
class TimeDomainEquationSystemProblemBuilder : public TimeDomainProblemBuilder,
                                               public EquationSystemProblemBuilderInterface
{
public:
  TimeDomainEquationSystemProblemBuilder()
    : TimeDomainProblemBuilder(nullptr),
      _problem{std::make_unique<TimeDomainEquationSystemProblem>()}
  {
  }

  ~TimeDomainEquationSystemProblemBuilder() override = default;

  void InitializeKernels() final;

  std::unique_ptr<TimeDomainEquationSystemProblem> ReturnProblem() { return std::move(_problem); }

protected:
  [[nodiscard]] TimeDomainEquationSystemProblem * GetProblem() const override
  {
    return _problem.get();
  }

  [[nodiscard]] TimeDependentEquationSystem * GetEquationSystem() const override
  {
    return GetProblem()->GetEquationSystem();
  }

private:
  std::unique_ptr<TimeDomainEquationSystemProblem> _problem{nullptr};
};

} // namespace hephaestus

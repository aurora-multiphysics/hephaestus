#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_builder.hpp"
#include "time_domain_equation_system_problem_operator.hpp"

namespace hephaestus
{
/// Time-depent problems with an equation system.
class TimeDomainEquationSystemProblem : public TimeDomainProblem,
                                        public EquationSystemProblemInterface
{
public:
  TimeDomainEquationSystemProblem() = default;
  ~TimeDomainEquationSystemProblem() override = default;

  [[nodiscard]] TimeDomainEquationSystemProblemOperator * GetOperator() const override
  {
    return static_cast<TimeDomainEquationSystemProblemOperator *>(TimeDomainProblem::GetOperator());
  }

  void SetOperator(std::unique_ptr<TimeDomainEquationSystemProblemOperator> problem_operator)
  {
    TimeDomainProblem::SetOperator(std::move(problem_operator));
  }

  void ConstructOperator() override
  {
    auto equation_system = std::make_unique<TimeDependentEquationSystem>();
    auto problem_operator = std::make_unique<TimeDomainEquationSystemProblemOperator>(
        *this, std::move(equation_system));

    SetOperator(std::move(problem_operator));
  }

  [[nodiscard]] TimeDependentEquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }
};

// Problem-builder for TimeDomainEquationSystemProblem.
class TimeDomainEquationSystemProblemBuilder : public TimeDomainProblemBuilder,
                                               public EquationSystemProblemBuilderInterface
{
public:
  /// NB: set "_problem" member variable in parent class.
  TimeDomainEquationSystemProblemBuilder()
    : TimeDomainProblemBuilder(std::make_unique<TimeDomainEquationSystemProblem>())
  {
  }

  ~TimeDomainEquationSystemProblemBuilder() override = default;

  /// NB: - note use of final. Ensure that the equation system is initialized.
  void InitializeKernels() final;

protected:
  /// NB: ensure @a GetProblem accessor is called in methods rather than using the "_problem" member variable.
  [[nodiscard]] TimeDomainEquationSystemProblem * GetProblem() const override
  {
    return static_cast<TimeDomainEquationSystemProblem *>(TimeDomainProblemBuilder::GetProblem());
  }

  [[nodiscard]] TimeDependentEquationSystem * GetEquationSystem() const override
  {
    return GetProblem()->GetEquationSystem();
  }
};

} // namespace hephaestus

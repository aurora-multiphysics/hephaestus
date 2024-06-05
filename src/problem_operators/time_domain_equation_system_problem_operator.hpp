#pragma once
#include "../common/pfem_extras.hpp"
#include "time_domain_problem_operator.hpp"
#include "problem_operator_base.hpp"
#include "equation_system_interface.hpp"

namespace hephaestus
{

/// Problem operator for time-dependent problems with an equation system.
class TimeDomainEquationSystemProblemOperator : public TimeDomainProblemOperator,
                                                public EquationSystemInterface
{
public:
  TimeDomainEquationSystemProblemOperator(const hephaestus::InputParameters &) = delete;
  TimeDomainEquationSystemProblemOperator(
      const hephaestus::InputParameters & params,
      std::unique_ptr<hephaestus::TimeDependentEquationSystem> equation_system)
    : TimeDomainProblemOperator(params), _equation_system{std::move(equation_system)}
  {
  }

  void Init() override;

  void Update() override;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;

  [[nodiscard]] hephaestus::TimeDependentEquationSystem * GetEquationSystem() const override
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added.");
    }

    return _equation_system.get();
  }

protected:
  void SetTrialVariableNames() override;
  void SetTrialVariables() override;

  void BuildEquationSystemOperator(double dt);

private:
  std::vector<mfem::ParGridFunction *> _trial_variable_time_derivatives;
  std::unique_ptr<hephaestus::TimeDependentEquationSystem> _equation_system{nullptr};
};

} // namespace hephaestus
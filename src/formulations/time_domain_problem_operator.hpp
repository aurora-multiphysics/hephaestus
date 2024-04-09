#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"
#include "problem_operator_interface.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

/// Problem operator for time-dependent problems with no equation system. The user will need to subclass this since the solve is not
/// implemented.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator,
                                  public ProblemOperatorInterface
{
public:
  TimeDomainProblemOperator(hephaestus::Problem & problem) : ProblemOperatorInterface(problem) {}
  ~TimeDomainProblemOperator() override = default;

  void SetGridFunctions() override {}
  void Init(mfem::Vector & X) override {}
  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override {}
};

/// Problem operator for time-dependent problems with an equation system.
class TimeDomainEquationSystemProblemOperator : public TimeDomainProblemOperator,
                                                public EquationSystemProblemOperatorInterface
{
public:
  TimeDomainEquationSystemProblemOperator(hephaestus::Problem &) = delete;
  TimeDomainEquationSystemProblemOperator(
      hephaestus::Problem & problem, std::unique_ptr<TimeDependentEquationSystem> equation_system)
    : TimeDomainProblemOperator(problem), _equation_system{std::move(equation_system)}
  {
  }

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;

  [[nodiscard]] TimeDependentEquationSystem * GetEquationSystem() const override
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added.");
    }

    return _equation_system.get();
  }

protected:
  void BuildEquationSystemOperator(double dt);

private:
  std::vector<mfem::ParGridFunction *> _trial_variable_time_derivatives;
  std::unique_ptr<TimeDependentEquationSystem> _equation_system{nullptr};
};

} // namespace hephaestus
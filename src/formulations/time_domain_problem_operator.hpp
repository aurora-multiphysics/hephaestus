#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"
#include "problem_operator_interface.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

/// Specifies output interfaces of a time-domain formulation. No equation-system! This class is not currently fully-implemented.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator,
                                  public ProblemOperatorInterface
{
public:
  TimeDomainProblemOperator(hephaestus::Problem & problem) : ProblemOperatorInterface(problem) {}

  void SetGridFunctions() override {}
  void Init(mfem::Vector & X) override {}

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override {}

protected:
  std::vector<mfem::ParGridFunction *> _trial_variable_time_derivatives;
};

/// TimeDomainProblemOperator with an equation system.
class TimeDomainEquationSystemProblemOperator : public TimeDomainProblemOperator,
                                                public EquationSystemProblemOperatorInterface
{
public:
  TimeDomainEquationSystemProblemOperator(hephaestus::Problem & problem)
    : TimeDomainProblemOperator(problem)
  {
  }

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;

  void SetEquationSystem(std::unique_ptr<TimeDependentEquationSystem> equation_system)
  {
    _equation_system.reset();
    _equation_system = std::move(equation_system);
  }

  [[nodiscard]] TimeDependentEquationSystem * GetEquationSystem() const override
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added to TimeDomainProblemOperator.");
    }

    return _equation_system.get();
  }

protected:
  void BuildEquationSystemOperator(double dt);

private:
  std::unique_ptr<TimeDependentEquationSystem> _equation_system{nullptr};
};

} // namespace hephaestus
#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"
#include "problem_operator_interface.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain formulation.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator,
                                  public EquationSystemProblemOperatorInterface
{
public:
  TimeDomainProblemOperator(hephaestus::Problem & problem)
    : EquationSystemProblemOperatorInterface(problem)
  {
  }

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
  void BuildEquationSystemOperator(double dt);

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
  std::unique_ptr<TimeDependentEquationSystem> _equation_system{nullptr};
  std::vector<mfem::ParGridFunction *> _trial_variable_time_derivatives;
};

} // namespace hephaestus
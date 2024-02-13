#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain formulation.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator
{
public:
  TimeDomainProblemOperator(hephaestus::Problem & problem) : _problem(problem) {}

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector & X);
  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
  virtual void BuildEquationSystemOperator(double dt);

  mfem::Array<int> _true_offsets, _block_true_offsets;
  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variable_time_derivatives, _trial_variables;

  mfem::BlockVector _true_x, _true_rhs;
  mfem::OperatorHandle _equation_system_operator;

  void SetEquationSystem(std::unique_ptr<TimeDependentEquationSystem> new_equation_system)
  {
    _equation_system.reset();
    _equation_system = std::move(new_equation_system);
  }

  inline EquationSystem * GetEquationSystem() const
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added to TimeDomainProblemOperator.");
    }

    return _equation_system.get();
  }

protected:
  hephaestus::Problem & _problem;
  std::unique_ptr<TimeDependentEquationSystem> _equation_system{nullptr};
};

} // namespace hephaestus
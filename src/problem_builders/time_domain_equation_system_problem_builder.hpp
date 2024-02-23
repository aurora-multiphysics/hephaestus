#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_operator.hpp"

namespace hephaestus
{
// Stores data required to describe a time domain formulation
class TimeDomainEquationSystemProblem : public hephaestus::EquationSystemProblem
{
public:
  friend class TimeDomainEquationSystemProblemBuilder;

  TimeDomainEquationSystemProblem() = default;
  ~TimeDomainEquationSystemProblem() override = default;

  [[nodiscard]] hephaestus::TimeDependentEquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }

  [[nodiscard]] hephaestus::TimeDomainEquationSystemProblemOperator * GetOperator() const override
  {
    if (!_td_operator)
    {
      MFEM_ABORT("No TimeDomainEquationSystemProblemOperator has been added to "
                 "TimeDomainEquationSystemProblem.");
    }

    return _td_operator.get();
  }

  void SetOperator(std::unique_ptr<TimeDomainEquationSystemProblemOperator> new_problem_operator)
  {
    _td_operator.reset();
    _td_operator = std::move(new_problem_operator);
  }

protected:
  std::unique_ptr<hephaestus::TimeDomainEquationSystemProblemOperator> _td_operator{nullptr};
};

// Builder class of a time-domain EM formulation.
class TimeDomainEquationSystemProblemBuilder : public EquationSystemProblemBuilder
{
public:
  TimeDomainEquationSystemProblemBuilder()
    : _problem(std::make_unique<hephaestus::TimeDomainEquationSystemProblem>())
  {
  }

  ~TimeDomainEquationSystemProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::TimeDomainEquationSystemProblem> ReturnProblem()
  {
    return std::move(_problem);
  }

  static std::vector<mfem::ParGridFunction *>
  RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                          hephaestus::GridFunctions & gridfunctions);

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override;

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override;

protected:
  std::unique_ptr<hephaestus::TimeDomainEquationSystemProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  hephaestus::TimeDomainEquationSystemProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus

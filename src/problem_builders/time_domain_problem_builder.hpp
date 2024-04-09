#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_operator.hpp"

namespace hephaestus
{

class TimeDomainProblem : public Problem
{
public:
  [[nodiscard]] TimeDomainProblemOperator * GetOperator() const override
  {
    if (!_problem_operator)
    {
      MFEM_ABORT("No operator has been added.");
    }

    return _problem_operator.get();
  }

  void SetOperator(std::unique_ptr<TimeDomainProblemOperator> problem_operator)
  {
    _problem_operator.reset();
    _problem_operator = std::move(problem_operator);
  }

  void ConstructOperator() override
  {
    _problem_operator.reset();
    _problem_operator = std::make_unique<TimeDomainProblemOperator>(*this);
  }

private:
  std::unique_ptr<TimeDomainProblemOperator> _problem_operator{nullptr};
};

class TimeDomainProblemBuilder : public ProblemBuilder
{
public:
  TimeDomainProblemBuilder() : _problem(std::make_unique<hephaestus::TimeDomainProblem>()) {}

  ~TimeDomainProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::TimeDomainProblem> ReturnProblem()
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
  std::unique_ptr<TimeDomainProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  TimeDomainProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus

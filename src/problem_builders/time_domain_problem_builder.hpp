#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_operator.hpp"

namespace hephaestus
{

// Stores data required to describe a time domain formulation
using TimeDomainProblem = ProblemTemplate<TimeDomainProblemOperator>;

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
  std::unique_ptr<hephaestus::TimeDomainProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  hephaestus::TimeDomainProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus

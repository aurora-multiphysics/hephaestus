#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus
{

class SteadyStateProblem : public hephaestus::Problem
{
public:
  std::unique_ptr<hephaestus::ProblemOperator> _ss_operator{nullptr};

  SteadyStateProblem() = default;
  ~SteadyStateProblem() override = default;

  hephaestus::EquationSystem * GetEquationSystem() override
  {
    return _ss_operator->GetEquationSystem();
  }

  hephaestus::ProblemOperator * GetOperator() override { return _ss_operator.get(); }
};

// Builder class of a frequency-domain problem.
class SteadyStateProblemBuilder : public hephaestus::ProblemBuilder
{
public:
  SteadyStateProblemBuilder() : _problem(std::make_unique<hephaestus::SteadyStateProblem>()) {}

  ~SteadyStateProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::SteadyStateProblem> ReturnProblem()
  {
    return std::move(_problem);
  }

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override {}

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void InitializeKernels() override;

  void ConstructEquationSystem() override {}

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

protected:
  std::unique_ptr<hephaestus::SteadyStateProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  hephaestus::SteadyStateProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus

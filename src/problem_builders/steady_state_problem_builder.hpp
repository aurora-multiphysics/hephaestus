#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus
{

/// Class for steady-state problems with no equation system.
class SteadyStateProblem : public Problem
{
public:
  SteadyStateProblem() = default;
  ~SteadyStateProblem() override = default;

  [[nodiscard]] ProblemOperator * GetOperator() const override
  {
    if (!_problem_operator)
    {
      MFEM_ABORT("No operator has been added.");
    }

    return _problem_operator.get();
  }

  void SetOperator(std::unique_ptr<ProblemOperator> problem_operator)
  {
    _problem_operator.reset();
    _problem_operator = std::move(problem_operator);
  }

  void ConstructOperator() override
  {
    _problem_operator.reset();
    _problem_operator = std::make_unique<ProblemOperator>(*this);
  }

private:
  std::unique_ptr<ProblemOperator> _problem_operator{nullptr};
};

class SteadyStateProblemBuilder : public ProblemBuilder
{
public:
  SteadyStateProblemBuilder(std::unique_ptr<SteadyStateProblem> problem)
    : _problem{std::move(problem)}
  {
  }
  SteadyStateProblemBuilder() : _problem(std::make_unique<hephaestus::SteadyStateProblem>()) {}

  ~SteadyStateProblemBuilder() override = default;

  std::unique_ptr<hephaestus::SteadyStateProblem> ReturnProblem() { return std::move(_problem); }

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override {}

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

protected:
  std::unique_ptr<SteadyStateProblem> _problem{nullptr};

  [[nodiscard]] SteadyStateProblem * GetProblem() const override { return _problem.get(); };

  mfem::ConstantCoefficient _one_coef{1.0};
};

} // namespace hephaestus

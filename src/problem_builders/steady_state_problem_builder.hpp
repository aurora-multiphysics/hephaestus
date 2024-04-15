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

  [[nodiscard]] hephaestus::ProblemOperator * GetOperator() const override
  {
    if (!_problem_operator)
    {
      MFEM_ABORT("No operator has been added.");
    }

    return _problem_operator.get();
  }

  void SetOperator(std::unique_ptr<hephaestus::ProblemOperator> problem_operator)
  {
    _problem_operator.reset();
    _problem_operator = std::move(problem_operator);
  }

  void ConstructOperator() override
  {
    _problem_operator.reset();
    _problem_operator = std::make_unique<hephaestus::ProblemOperator>(*this);
  }

private:
  std::unique_ptr<hephaestus::ProblemOperator> _problem_operator{nullptr};
};

class SteadyStateProblemBuilder : public ProblemBuilder
{
public:
  /// NB: constructor called in derived classes. The problem must be a subclass of SteadyStateProblem.
  SteadyStateProblemBuilder(std::unique_ptr<hephaestus::SteadyStateProblem> problem)
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
  /// NB: it is extremely important that this accessor is used in all member variables since
  /// derived classes will override this method. This allows us to reuse the methods defined
  /// here without having to override them in derived classes. Calling "_problem" directly
  /// will get you into trouble (likely a segfault since it will be NULL if this is a parent class.)
  [[nodiscard]] hephaestus::SteadyStateProblem * GetProblem() const override
  {
    return _problem.get();
  };

private:
  std::unique_ptr<hephaestus::SteadyStateProblem> _problem{nullptr};
};

} // namespace hephaestus

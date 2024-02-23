#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus
{

/// Class for steady-state problems with an equation system.
class SteadyStateEquationSystemProblem : public EquationSystemProblem
{
public:
  friend class SteadyStateEquationSystemProblemBuilder;

  SteadyStateEquationSystemProblem() = default;
  ~SteadyStateEquationSystemProblem() override = default;

  [[nodiscard]] hephaestus::EquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }

  [[nodiscard]] hephaestus::EquationSystemProblemOperator * GetOperator() const override
  {
    if (!_ss_operator)
    {
      MFEM_ABORT("No ProblemOperator has been added to SteadyStateProblem.");
    }

    return _ss_operator.get();
  }

  void SetOperator(std::unique_ptr<EquationSystemProblemOperator> new_problem_operator)
  {
    _ss_operator.reset();
    _ss_operator = std::move(new_problem_operator);
  }

private:
  std::unique_ptr<hephaestus::EquationSystemProblemOperator> _ss_operator{nullptr};
};

// Builder class of a frequency-domain problem.
class SteadyStateEquationSystemProblemBuilder : public EquationSystemProblemBuilder
{
public:
  SteadyStateEquationSystemProblemBuilder()
    : _problem(std::make_unique<hephaestus::SteadyStateEquationSystemProblem>())
  {
  }

  ~SteadyStateEquationSystemProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::SteadyStateEquationSystemProblem> ReturnProblem()
  {
    return std::move(_problem);
  }

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override {}

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

protected:
  std::unique_ptr<hephaestus::SteadyStateEquationSystemProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  hephaestus::SteadyStateEquationSystemProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus

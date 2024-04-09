#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

/// Class for steady-state problems with an equation system.
class SteadyStateEquationSystemProblem : public SteadyStateProblem,
                                         public EquationSystemProblemInterface
{
public:
  [[nodiscard]] EquationSystemProblemOperator * GetOperator() const override
  {
    if (!_problem_operator)
    {
      MFEM_ABORT("No operator has been added.");
    }

    return _problem_operator.get();
  }

  void SetOperator(std::unique_ptr<EquationSystemProblemOperator> problem_operator)
  {
    _problem_operator.reset();
    _problem_operator = std::move(problem_operator);
  }

  void ConstructOperator() override
  {
    auto equation_system = std::make_unique<EquationSystem>();

    _problem_operator.reset();
    _problem_operator =
        std::make_unique<EquationSystemProblemOperator>(*this, std::move(equation_system));
  }

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    return GetOperator()->GetEquationSystem();
  }

private:
  std::unique_ptr<EquationSystemProblemOperator> _problem_operator{nullptr};
};

class SteadyStateEquationSystemProblemBuilder : public SteadyStateProblemBuilder,
                                                public EquationSystemProblemBuilderInterface
{
public:
  SteadyStateEquationSystemProblemBuilder()
    : SteadyStateProblemBuilder(nullptr),
      _problem{std::make_unique<SteadyStateEquationSystemProblem>()}
  {
  }

  void InitializeKernels() final
  {
    ProblemBuilder::InitializeKernels();

    GetEquationSystem()->Init(GetProblem()->_gridfunctions,
                              GetProblem()->_fespaces,
                              GetProblem()->_bc_map,
                              GetProblem()->_coefficients);
  }

  std::unique_ptr<SteadyStateEquationSystemProblem> ReturnProblem() { return std::move(_problem); }

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override {}

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

protected:
  std::unique_ptr<SteadyStateEquationSystemProblem> _problem{nullptr};

  [[nodiscard]] SteadyStateEquationSystemProblem * GetProblem() const override
  {
    return _problem.get();
  }

  [[nodiscard]] EquationSystem * GetEquationSystem() const override
  {
    return GetProblem()->GetEquationSystem();
  }
};

} // namespace hephaestus

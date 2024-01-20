#pragma once
#include "equation_system_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus
{

class SteadyStateProblem : public hephaestus::Problem
{
public:
  std::unique_ptr<hephaestus::EquationSystem> eq_sys{nullptr};
  std::unique_ptr<hephaestus::EquationSystemOperator> eq_sys_operator{nullptr};

  SteadyStateProblem() = default;
  ~SteadyStateProblem() override = default;

  hephaestus::EquationSystem * GetEquationSystem() override { return eq_sys.get(); };
  hephaestus::EquationSystemOperator * GetOperator() override { return eq_sys_operator.get(); };
};

// Builder class of a frequency-domain problem.
class SteadyStateProblemBuilder : public hephaestus::ProblemBuilder
{
public:
  SteadyStateProblemBuilder() : problem(std::make_unique<hephaestus::SteadyStateProblem>()){};

  ~SteadyStateProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::SteadyStateProblem> ReturnProblem()
  {
    return std::move(problem);
  };

  void RegisterFESpaces() override{};

  void RegisterGridFunctions() override{};

  void RegisterAuxSolvers() override{};

  void RegisterCoefficients() override{};

  void InitializeKernels() override;

  void ConstructEquationSystem() override{};

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructSolver() override{};

protected:
  std::unique_ptr<hephaestus::SteadyStateProblem> problem{nullptr};
  mfem::ConstantCoefficient oneCoef{1.0};

  hephaestus::SteadyStateProblem * GetProblem() override { return problem.get(); };
};

} // namespace hephaestus

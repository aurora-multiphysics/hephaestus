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
  ~SteadyStateProblem() override {}

  virtual hephaestus::EquationSystem * GetEquationSystem() { return eq_sys.get(); };
  virtual hephaestus::EquationSystemOperator * GetOperator() { return eq_sys_operator.get(); };
};

// Builder class of a frequency-domain problem.
class SteadyStateProblemBuilder : public hephaestus::ProblemBuilder
{
public:
  SteadyStateProblemBuilder() : problem(std::make_unique<hephaestus::SteadyStateProblem>()){};

  ~SteadyStateProblemBuilder() override {}

  virtual std::unique_ptr<hephaestus::SteadyStateProblem> ReturnProblem()
  {
    return std::move(problem);
  };

  virtual void RegisterFESpaces() override{};

  virtual void RegisterGridFunctions() override{};

  virtual void RegisterAuxSolvers() override{};

  virtual void RegisterCoefficients() override{};

  virtual void InitializeKernels() override;

  virtual void ConstructEquationSystem() override{};

  virtual void ConstructOperator() override;

  virtual void ConstructState() override;

  virtual void ConstructSolver() override{};

protected:
  std::unique_ptr<hephaestus::SteadyStateProblem> problem{nullptr};
  mfem::ConstantCoefficient oneCoef{1.0};

  virtual hephaestus::SteadyStateProblem * GetProblem() override { return problem.get(); };
};

} // namespace hephaestus

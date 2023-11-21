#pragma once
#include "problem_builder_base.hpp"
#include "problem_operator.hpp"
namespace hephaestus {

class SteadyStateProblem : public hephaestus::Problem {
public:
  std::unique_ptr<hephaestus::EquationSystem> eq_sys;
  std::unique_ptr<hephaestus::ProblemOperator> eq_sys_operator;

  SteadyStateProblem() = default;

  virtual hephaestus::EquationSystem *GetEquationSystem() {
    return eq_sys.get();
  };
  virtual hephaestus::ProblemOperator *GetOperator() {
    return eq_sys_operator.get();
  };
};

// Builder class of a frequency-domain problem.
class SteadyStateProblemBuilder : public hephaestus::ProblemBuilder {
protected:
  std::unique_ptr<hephaestus::SteadyStateProblem> problem;
  mfem::ConstantCoefficient oneCoef{1.0};

  virtual hephaestus::SteadyStateProblem *GetProblem() override {
    return this->problem.get();
  };

public:
  SteadyStateProblemBuilder()
      : problem(std::make_unique<hephaestus::SteadyStateProblem>()){};

  virtual std::unique_ptr<hephaestus::SteadyStateProblem> ReturnProblem() {
    return std::move(this->problem);
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
};

} // namespace hephaestus

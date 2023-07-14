#pragma once
#include "frequency_domain_equation_system_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus {

/** Class that store objects associated with a set of PDEs to be solved in the
 * frequency domain. */
class FrequencyDomainProblem : public hephaestus::Problem {
public:
  std::unique_ptr<hephaestus::EquationSystem> eq_sys;
  std::unique_ptr<hephaestus::FrequencyDomainEquationSystemOperator>
      fd_operator;

  FrequencyDomainProblem() = default;

  virtual hephaestus::EquationSystem *GetEquationSystem() {
    return eq_sys.get();
  };
  virtual hephaestus::FrequencyDomainEquationSystemOperator *GetOperator() {
    return fd_operator.get();
  };
};

// Builder class of a frequency-domain problem.
class FrequencyDomainProblemBuilder : public hephaestus::ProblemBuilder {
protected:
  std::unique_ptr<hephaestus::FrequencyDomainProblem> problem;
  mfem::ConstantCoefficient oneCoef{1.0};
  mfem::ConstantCoefficient *freqCoef;

  virtual hephaestus::FrequencyDomainProblem *GetProblem() override {
    return this->problem.get();
  };

public:
  FrequencyDomainProblemBuilder()
      : problem(std::make_unique<hephaestus::FrequencyDomainProblem>()){};

  virtual std::unique_ptr<hephaestus::FrequencyDomainProblem> ReturnProblem() {
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

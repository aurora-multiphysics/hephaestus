#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_operator.hpp"

namespace hephaestus {

// Stores data required to describe a time domain formulation
class TimeDomainProblem : public hephaestus::Problem {
public:
  std::unique_ptr<hephaestus::TimeDependentEquationSystem> td_equation_system;
  std::unique_ptr<hephaestus::TimeDomainProblemOperator> td_operator;

  TimeDomainProblem() = default;

  virtual hephaestus::TimeDependentEquationSystem *GetEquationSystem() {
    return td_equation_system.get();
  };
  virtual hephaestus::TimeDomainProblemOperator *GetOperator() {
    return td_operator.get();
  };
};

// Builder class of a time-domain EM formulation.
class TimeDomainProblemBuilder : public hephaestus::ProblemBuilder {
protected:
  std::unique_ptr<hephaestus::TimeDomainProblem> problem;
  mfem::ConstantCoefficient oneCoef{1.0};

  virtual hephaestus::TimeDomainProblem *GetProblem() override {
    return this->problem.get();
  };

public:
  TimeDomainProblemBuilder()
      : problem(std::make_unique<hephaestus::TimeDomainProblem>()){};

  virtual std::unique_ptr<hephaestus::TimeDomainProblem> ReturnProblem() {
    return std::move(this->problem);
  };

  static std::vector<mfem::ParGridFunction *>
  RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                          hephaestus::GridFunctions &gridfunctions);

  virtual void RegisterFESpaces() override{};

  virtual void RegisterGridFunctions() override;

  virtual void RegisterAuxSolvers() override{};

  virtual void RegisterCoefficients() override{};

  virtual void ConstructEquationSystem() override;

  virtual void ConstructJacobianPreconditioner() override;

  virtual void ConstructJacobianSolver() override;

  virtual void ConstructNonlinearSolver() override;

  virtual void InitializeKernels() override;

  virtual void ConstructOperator() override;

  virtual void ConstructState() override;

  virtual void ConstructTimestepper() override;
};

} // namespace hephaestus

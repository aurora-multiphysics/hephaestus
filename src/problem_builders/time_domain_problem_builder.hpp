#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_equation_system_operator.hpp"

namespace hephaestus
{

// Stores data required to describe a time domain formulation
class TimeDomainProblem : public hephaestus::Problem
{
public:
  std::unique_ptr<hephaestus::TimeDependentEquationSystem> td_equation_system;
  std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator> td_operator;

  TimeDomainProblem() = default;
  ~TimeDomainProblem() override;

  hephaestus::TimeDependentEquationSystem * GetEquationSystem() override
  {
    return td_equation_system.get();
  };
  hephaestus::TimeDomainEquationSystemOperator * GetOperator() override
  {
    return td_operator.get();
  };
};

// Builder class of a time-domain EM formulation.
class TimeDomainProblemBuilder : public hephaestus::ProblemBuilder
{
public:
  TimeDomainProblemBuilder() : problem(std::make_unique<hephaestus::TimeDomainProblem>()){};

  ~TimeDomainProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::TimeDomainProblem> ReturnProblem()
  {
    return std::move(problem);
  };

  static std::vector<mfem::ParGridFunction *>
  RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                          hephaestus::GridFunctions & gridfunctions);

  void RegisterFESpaces() override{};

  void RegisterGridFunctions() override;

  void RegisterAuxSolvers() override{};

  void RegisterCoefficients() override{};

  void ConstructEquationSystem() override;

  void InitializeKernels() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructSolver() override;

protected:
  std::unique_ptr<hephaestus::TimeDomainProblem> problem{nullptr};
  mfem::ConstantCoefficient oneCoef{1.0};

  hephaestus::TimeDomainProblem * GetProblem() override { return problem.get(); };
};

} // namespace hephaestus

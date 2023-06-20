#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_equation_system_operator.hpp"

namespace hephaestus {

// Stores data required to describe a time domain formulation
class TimeDomainProblem : public hephaestus::Problem {
public:
  std::unique_ptr<hephaestus::TimeDependentEquationSystem> td_equation_system;
  std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator> td_operator;

  TimeDomainProblem() = default;

  virtual hephaestus::TimeDependentEquationSystem *GetEquationSystem() {
    return td_equation_system.get();
  };
  virtual hephaestus::TimeDomainEquationSystemOperator *GetOperator() {
    return td_operator.get();
  };
};

// Builder class of a time-domain EM formulation.
class TimeDomainProblemBuilder : public hephaestus::ProblemBuilder {
protected:
  std::unique_ptr<hephaestus::TimeDomainProblem> problem;
  virtual hephaestus::TimeDomainProblem *GetProblem() override {
    return this->problem.get();
  };

public:
  TimeDomainProblemBuilder()
      : problem(std::make_unique<hephaestus::TimeDomainProblem>()){};
  mfem::ConstantCoefficient oneCoef{1.0};
  virtual std::unique_ptr<hephaestus::TimeDomainProblem> ReturnProblem() {
    return std::move(this->problem);
  };

  virtual std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator>
  CreateTimeDomainEquationSystemOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) const;

  static std::vector<mfem::ParGridFunction *> RegisterTimeDerivatives(
      std::vector<std::string> gridfunction_names,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &gridfunctions);

  virtual void RegisterFESpaces() override{};

  virtual void RegisterGridFunctions() override {
    std::vector<std::string> gridfunction_names;
    for (auto const &[name, gf] : this->problem->gridfunctions) {
      gridfunction_names.push_back(name);
    }
    this->RegisterTimeDerivatives(gridfunction_names,
                                  this->problem->gridfunctions);
  };

  virtual void RegisterAuxSolvers() override{};

  virtual void RegisterCoefficients() override{};

  virtual void ConstructEquationSystem() override {
    hephaestus::InputParameters params;
    this->problem->td_equation_system =
        std::make_unique<hephaestus::TimeDependentEquationSystem>(params);
  };

  virtual void InitializeKernels() override {
    this->problem->td_equation_system->Init(
        this->problem->gridfunctions, this->problem->fespaces,
        this->problem->bc_map, this->problem->domain_properties);

    this->problem->preprocessors.Init(this->problem->gridfunctions,
                                      this->problem->domain_properties);
    this->problem->sources.Init(this->problem->gridfunctions,
                                this->problem->fespaces, this->problem->bc_map,
                                this->problem->domain_properties);
  };

  virtual void ConstructOperator() override {
    this->problem->td_operator = this->CreateTimeDomainEquationSystemOperator(
        *(this->problem->pmesh), this->problem->fespaces,
        this->problem->gridfunctions, this->problem->bc_map,
        this->problem->domain_properties, this->problem->sources,
        this->problem->solver_options);
    this->problem->td_operator->SetEquationSystem(
        this->problem->td_equation_system.get());
    this->problem->td_operator->SetVariables();
  };

  virtual void ConstructState() override {
    this->problem->F = new mfem::BlockVector(
        this->problem->td_operator->true_offsets); // Vector of dofs
    this->problem->td_operator->Init(
        *(this->problem->F)); // Set up initial conditions
    this->problem->td_operator->SetTime(0.0);
  };

  virtual void ConstructSolver() override {
    this->problem->ode_solver = new mfem::BackwardEulerSolver;
    this->problem->ode_solver->Init(*(this->problem->td_operator));
  };
};

} // namespace hephaestus

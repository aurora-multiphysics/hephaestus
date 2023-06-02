#pragma once
#include "problem_builder_base.hpp"

namespace hephaestus {

// Stores data required to describe a time domain formulation
class TransientProblem : public hephaestus::Problem {
public:
  hephaestus::TransientFormulation *formulation;
  std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator> td_operator;
  std::unique_ptr<hephaestus::TimeDependentEquationSystem> td_equation_system;

  TransientProblem() = default;
  explicit TransientProblem(const hephaestus::InputParameters &params);

  virtual hephaestus::TransientFormulation *GetFormulation() {
    return formulation;
  };
  virtual hephaestus::TimeDependentEquationSystem *GetEquationSystem() {
    return td_equation_system.get();
  };
  virtual hephaestus::TimeDomainEquationSystemOperator *GetOperator() {
    return td_operator.get();
  };
};

// Builder class of a time-domain EM formulation.
class TransientProblemBuilder : public hephaestus::ProblemBuilder {
private:
  std::unique_ptr<hephaestus::TransientProblem> problem;
  virtual hephaestus::TransientProblem *GetProblem() override {
    return this->problem.get();
  };

public:
  TransientProblemBuilder()
      : problem(std::make_unique<hephaestus::TransientProblem>()){};
  TransientProblemBuilder(const hephaestus::InputParameters &params)
      : problem(std::make_unique<hephaestus::TransientProblem>(params)){};

  virtual std::unique_ptr<hephaestus::TransientProblem> ReturnProblem() {
    return std::move(this->problem);
  };

  virtual void SetFormulation(hephaestus::TransientFormulation *formulation) {
    this->problem->formulation = formulation;
    ;
  }

  virtual void RegisterFESpaces() override {
    this->problem->fespaces.Init(*(this->problem->pmesh));
  };

  virtual void RegisterGridFunctions() override {
    this->problem->gridfunctions.Init(*(this->problem->pmesh),
                                      this->problem->fespaces);
    std::vector<std::string> gridfunction_names;

    for (auto const &[name, gf] : this->problem->gridfunctions) {
      gridfunction_names.push_back(name);
    }
    this->problem->formulation->RegisterTimeDerivatives(
        gridfunction_names, this->problem->gridfunctions);

    this->problem->formulation->RegisterMissingVariables(
        *(this->problem->pmesh), this->problem->fespaces,
        this->problem->gridfunctions);
  };

  virtual void RegisterAuxKernels() override {
    this->problem->formulation->RegisterAuxKernels(this->problem->gridfunctions,
                                                   this->problem->auxkernels);
  };

  virtual void RegisterCoefficients() override {
    this->problem->formulation->RegisterCoefficients(
        this->problem->domain_properties);
  };

  virtual void ConstructEquationSystem() override {
    this->problem->td_equation_system =
        this->problem->formulation->CreateTimeDependentEquationSystem();
  };

  virtual void InitializeKernels() override {
    this->problem->td_equation_system->Init(
        this->problem->gridfunctions, this->problem->fespaces,
        this->problem->bc_map, this->problem->domain_properties);

    this->problem->auxkernels.Init(this->problem->gridfunctions,
                                   this->problem->domain_properties);
    this->problem->sources.Init(this->problem->gridfunctions,
                                this->problem->fespaces, this->problem->bc_map,
                                this->problem->domain_properties);
  };

  virtual void ConstructOperator() override {
    this->problem->td_operator =
        this->problem->formulation->CreateTimeDomainEquationSystemOperator(
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

  virtual void InitializePostprocessors() override {
    this->problem->postprocessors.Init(this->problem->gridfunctions,
                                       this->problem->domain_properties);
    this->problem->auxkernels.Solve(0.0);
  };
};

} // namespace hephaestus

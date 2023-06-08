#pragma once
#include "problem_builder_base.hpp"

namespace hephaestus {

class FrequencyDomainProblem : public hephaestus::Problem {
public:
  hephaestus::SteadyFormulation *formulation;
  std::unique_ptr<hephaestus::EquationSystem> eq_sys;
  std::unique_ptr<hephaestus::FrequencyDomainOperator> fd_operator;

  FrequencyDomainProblem() = default;
  explicit FrequencyDomainProblem(const hephaestus::InputParameters &params);

  virtual hephaestus::SteadyFormulation *GetFormulation() {
    return formulation;
  };
  virtual hephaestus::EquationSystem *GetEquationSystem() {
    return eq_sys.get();
  };
  virtual hephaestus::FrequencyDomainOperator *GetOperator() {
    return fd_operator.get();
  };
};

// Builder class of a frequency-domain problem.
class FrequencyDomainProblemBuilder : public hephaestus::ProblemBuilder {
private:
  std::unique_ptr<hephaestus::FrequencyDomainProblem> problem;

  virtual hephaestus::FrequencyDomainProblem *GetProblem() override {
    return this->problem.get();
  };

public:
  FrequencyDomainProblemBuilder()
      : problem(std::make_unique<hephaestus::FrequencyDomainProblem>()){};
  FrequencyDomainProblemBuilder(const hephaestus::InputParameters &params)
      : hephaestus::ProblemBuilder(params),
        problem(std::make_unique<hephaestus::FrequencyDomainProblem>(params)){};

  virtual std::unique_ptr<hephaestus::FrequencyDomainProblem> ReturnProblem() {
    return std::move(this->problem);
  };

  virtual void SetFormulation(hephaestus::SteadyFormulation *formulation) {
    this->problem->formulation = formulation;
  }

  virtual void RegisterFESpaces() override{
      // this->problem->fespaces.Init(*(this->problem->pmesh));
  };

  virtual void RegisterGridFunctions() override {
    // this->problem->gridfunctions.Init(*(this->problem->pmesh),
    //                                   this->problem->fespaces);
    this->problem->formulation->RegisterMissingVariables(
        *(this->problem->pmesh), this->problem->fespaces,
        this->problem->gridfunctions);
  };

  virtual void RegisterAuxSolvers() override {
    this->problem->formulation->RegisterAuxSolvers(
        this->problem->gridfunctions, this->problem->postprocessors);
  };

  virtual void RegisterCoefficients() override {
    this->problem->formulation->RegisterCoefficients(
        this->problem->domain_properties);
  };

  virtual void InitializeKernels() override {
    this->problem->preprocessors.Init(this->problem->gridfunctions,
                                      this->problem->domain_properties);
    this->problem->sources.Init(this->problem->gridfunctions,
                                this->problem->fespaces, this->problem->bc_map,
                                this->problem->domain_properties);
  };

  virtual void ConstructEquationSystem() override{};

  virtual void ConstructOperator() override {
    this->problem->fd_operator =
        this->problem->formulation->CreateFrequencyDomainOperator(
            *(this->problem->pmesh), this->problem->fespaces,
            this->problem->gridfunctions, this->problem->bc_map,
            this->problem->domain_properties, this->problem->sources,
            this->problem->solver_options);
    this->problem->fd_operator->SetVariables();
  };

  virtual void ConstructState() override {
    this->problem->F = new mfem::BlockVector(
        this->problem->fd_operator->true_offsets); // Vector of dofs
    this->problem->fd_operator->Init(
        *(this->problem->F)); // Set up initial conditions
  };

  virtual void ConstructSolver() override{};

  virtual void InitializePostprocessors() override {
    this->problem->postprocessors.Init(this->problem->gridfunctions,
                                       this->problem->domain_properties);
  };
};

} // namespace hephaestus
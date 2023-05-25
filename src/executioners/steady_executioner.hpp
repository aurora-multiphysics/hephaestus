#pragma once
#include "executioner_base.hpp"

namespace hephaestus {

class FrequencyDomainProblem : public hephaestus::Problem {
public:
  hephaestus::SteadyFormulation *formulation;
  std::unique_ptr<hephaestus::FrequencyDomainOperator> fd_operator;
  hephaestus::EquationSystem *eq_sys;

  FrequencyDomainProblem() = default;
  explicit FrequencyDomainProblem(const hephaestus::InputParameters &params);
};

// Builder class of a time-domain EM formulation.
class FrequencyDomainProblemBuilder {
private:
  std::unique_ptr<hephaestus::FrequencyDomainProblem> problem;

public:
  FrequencyDomainProblemBuilder();
  FrequencyDomainProblemBuilder(const hephaestus::InputParameters &params)
      : problem(std::make_unique<hephaestus::FrequencyDomainProblem>(params)){};

  virtual std::unique_ptr<hephaestus::FrequencyDomainProblem> GetProblem() {
    return std::move(this->problem);
  };

  virtual void RegisterFESpaces() {
    this->problem->fespaces.Init(this->problem->pmesh);
  };

  virtual void RegisterGridFunctions() {
    this->problem->gridfunctions.Init(this->problem->pmesh,
                                      this->problem->fespaces);
    this->problem->formulation->RegisterMissingVariables(
        this->problem->pmesh, this->problem->fespaces,
        this->problem->gridfunctions);
  };

  virtual void RegisterAuxKernels() {
    this->problem->formulation->RegisterAuxKernels(this->problem->gridfunctions,
                                                   this->problem->auxkernels);
  };

  virtual void RegisterCoefficients() {
    this->problem->formulation->RegisterCoefficients(
        this->problem->domain_properties);
  };

  virtual void InitializeKernels() {
    this->problem->auxkernels.Init(this->problem->gridfunctions,
                                   this->problem->domain_properties);
    this->problem->sources.Init(this->problem->gridfunctions,
                                this->problem->fespaces, this->problem->bc_map,
                                this->problem->domain_properties);
  };

  virtual void ConstructOperator() {
    this->problem->fd_operator =
        this->problem->formulation->CreateFrequencyDomainOperator(
            this->problem->pmesh, this->problem->fespaces,
            this->problem->gridfunctions, this->problem->bc_map,
            this->problem->domain_properties, this->problem->sources,
            this->problem->solver_options);
    this->problem->fd_operator->SetVariables();
  };

  virtual void ConstructState() {
    this->problem->F = new mfem::BlockVector(
        this->problem->fd_operator->true_offsets); // Vector of dofs
    this->problem->fd_operator->Init(
        *(this->problem->F)); // Set up initial conditions
  };

  virtual void InitializePostprocessors() {
    this->problem->postprocessors.Init(this->problem->gridfunctions,
                                       this->problem->domain_properties);
    this->problem->auxkernels.Solve();
  };
};

class SteadyExecutioner : public Executioner {
private:
  std::unique_ptr<hephaestus::FrequencyDomainProblemBuilder> problem_builder;

public:
  SteadyExecutioner() = default;
  explicit SteadyExecutioner(const hephaestus::InputParameters &params);

  void Init() override;

  void Solve() const override;

  void Execute() const override;
  std::unique_ptr<hephaestus::FrequencyDomainProblem> problem;
};

} // namespace hephaestus

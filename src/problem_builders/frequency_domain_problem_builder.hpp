#pragma once
#include "frequency_domain_equation_system_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus {

class FrequencyDomainProblem : public hephaestus::Problem {
public:
  std::unique_ptr<hephaestus::EquationSystem> eq_sys;
  std::unique_ptr<hephaestus::FrequencyDomainEquationSystemOperator>
      fd_operator;

  FrequencyDomainProblem() = default;
  explicit FrequencyDomainProblem(const hephaestus::InputParameters &params);

  virtual hephaestus::EquationSystem *GetEquationSystem() {
    return eq_sys.get();
  };
  virtual hephaestus::FrequencyDomainEquationSystemOperator *GetOperator() {
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
  mfem::ConstantCoefficient oneCoef{1.0};
  mfem::ConstantCoefficient *freqCoef;
  virtual std::unique_ptr<hephaestus::FrequencyDomainProblem> ReturnProblem() {
    return std::move(this->problem);
  };

  virtual std::unique_ptr<hephaestus::FrequencyDomainEquationSystemOperator>
  CreateFrequencyDomainEquationSystemOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) const {
    return std::make_unique<hephaestus::FrequencyDomainEquationSystemOperator>(
        pmesh, fespaces, variables, bc_map, domain_properties, sources,
        solver_options);
  };

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables){};

  virtual void
  RegisterAuxSolvers(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxSolvers &auxsolvers){};

  virtual void
  RegisterCoefficients(hephaestus::DomainProperties &domain_properties){};

  virtual void RegisterFESpaces() override{};

  virtual void RegisterGridFunctions() override {
    this->RegisterMissingVariables(*(this->problem->pmesh),
                                   this->problem->fespaces,
                                   this->problem->gridfunctions);
  };

  virtual void RegisterAuxSolvers() override {
    this->RegisterAuxSolvers(this->problem->gridfunctions,
                             this->problem->postprocessors);
  };

  virtual void RegisterCoefficients() override {
    this->RegisterCoefficients(this->problem->domain_properties);
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
        this->CreateFrequencyDomainEquationSystemOperator(
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

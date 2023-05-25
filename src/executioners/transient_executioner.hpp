#pragma once
#include "executioner_base.hpp"

namespace hephaestus {

// Stores data required to describe a time domain formulation
class TransientProblem : public hephaestus::Problem {
public:
  hephaestus::TransientFormulation *formulation;
  std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator> td_operator;
  std::unique_ptr<hephaestus::TimeDependentEquationSystem> td_equation_system;

  TransientProblem() = default;
  explicit TransientProblem(const hephaestus::InputParameters &params);
};

// Builder class of a time-domain EM formulation.
class TransientProblemBuilder {
private:
  std::unique_ptr<hephaestus::TransientProblem> problem;

public:
  TransientProblemBuilder();
  TransientProblemBuilder(const hephaestus::InputParameters &params)
      : problem(std::make_unique<hephaestus::TransientProblem>(params)){};

  virtual std::unique_ptr<hephaestus::TransientProblem> GetProblem() {
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

  virtual void ConstructEquationSystem() {
    this->problem->td_equation_system->Init(
        this->problem->gridfunctions, this->problem->fespaces,
        this->problem->bc_map, this->problem->domain_properties);
  };

  virtual void InitializeKernels() {
    this->problem->auxkernels.Init(this->problem->gridfunctions,
                                   this->problem->domain_properties);
    this->problem->sources.Init(this->problem->gridfunctions,
                                this->problem->fespaces, this->problem->bc_map,
                                this->problem->domain_properties);
  };

  virtual void ConstructOperator() {
    this->problem->td_operator =
        this->problem->formulation->CreateTimeDomainEquationSystemOperator(
            this->problem->pmesh, this->problem->fespaces,
            this->problem->gridfunctions, this->problem->bc_map,
            this->problem->domain_properties, this->problem->sources,
            this->problem->solver_options);
    this->problem->td_operator->SetEquationSystem(
        this->problem->td_equation_system.get());
    this->problem->td_operator->SetVariables();
  };

  virtual void ConstructState() {
    this->problem->F = new mfem::BlockVector(
        this->problem->td_operator->true_offsets); // Vector of dofs
    this->problem->td_operator->Init(
        *(this->problem->F)); // Set up initial conditions
    this->problem->td_operator->SetTime(0.0);
  };

  virtual void ConstructSolver() {
    this->problem->ode_solver = new mfem::BackwardEulerSolver;
    this->problem->ode_solver->Init(*(this->problem->td_operator));
  };

  virtual void InitializePostprocessors() {
    this->problem->postprocessors.Init(this->problem->gridfunctions,
                                       this->problem->domain_properties);
    this->problem->auxkernels.Solve(0.0);
  };
};

class TransientExecutioner : public Executioner {
private:
  double t_initial;       // Start time
  double t_final;         // End time
  mutable double t;       // Current time
  mutable int it;         // Time index
  int vis_steps;          // Number of cyces between each output update
  mutable bool last_step; // Flag to check if current step is final
  std::unique_ptr<hephaestus::TransientProblemBuilder> problem_builder;

public:
  mutable double t_step; // Time step

  TransientExecutioner() = default;
  explicit TransientExecutioner(const hephaestus::InputParameters &params);

  void Init() override;

  void Step(double dt, int it) const;

  void Solve() const override;

  void Execute() const override;
  std::unique_ptr<hephaestus::TransientProblem> problem;
};

} // namespace hephaestus

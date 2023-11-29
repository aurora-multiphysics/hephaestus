#include "steady_state_problem_builder.hpp"

namespace hephaestus {

void SteadyStateProblemBuilder::InitializeKernels() {
  this->problem->preprocessors.Init(this->problem->gridfunctions,
                                    this->problem->coefficients);
  this->problem->sources.Init(this->problem->gridfunctions,
                              this->problem->fespaces, this->problem->bc_map,
                              this->problem->coefficients);
}
void SteadyStateProblemBuilder::ConstructOperator() {
  this->problem->eq_sys_operator =
      std::make_unique<hephaestus::EquationSystemOperator>(
          *(this->problem->pmesh), this->problem->fespaces,
          this->problem->gridfunctions, this->problem->bc_map,
          this->problem->coefficients, this->problem->sources,
          this->problem->solver_options);
  this->problem->eq_sys_operator->SetGridFunctions();
}

void SteadyStateProblemBuilder::ConstructState() {
  this->problem->F = new mfem::BlockVector(
      this->problem->eq_sys_operator->true_offsets); // Vector of dofs
  this->problem->eq_sys_operator->Init(
      *(this->problem->F)); // Set up initial conditions
}

void SteadyStateProblemBuilder::ConstructEquationSystem() {
  hephaestus::InputParameters params;
  this->problem->eq_sys = 
      std::make_unique<hephaestus::EquationSystem>(params);
}
} // namespace hephaestus

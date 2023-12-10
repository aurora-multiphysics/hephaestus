#include "steady_state_problem_builder.hpp"

namespace hephaestus {

void SteadyStateProblemBuilder::InitializeKernels() {
  problem->preprocessors.Init(problem->gridfunctions, problem->coefficients);
  problem->sources.Init(problem->gridfunctions, problem->fespaces,
                        problem->bc_map, problem->coefficients);
}

void SteadyStateProblemBuilder::ConstructOperator() {
  problem->ss_operator =
      std::make_unique<hephaestus::ProblemOperator>(*problem);
  problem->ss_operator->SetGridFunctions();
}

void SteadyStateProblemBuilder::ConstructState() {
  problem->F = new mfem::BlockVector(
      problem->ss_operator->true_offsets);   // Vector of dofs
  problem->ss_operator->Init(*(problem->F)); // Set up initial conditions
}
} // namespace hephaestus

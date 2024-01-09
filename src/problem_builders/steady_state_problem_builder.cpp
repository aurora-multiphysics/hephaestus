#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::InitializeKernels()
{
  problem->preprocessors.Init(problem->gridfunctions, problem->coefficients);
  problem->sources.Init(
      problem->gridfunctions, problem->fespaces, problem->bc_map, problem->coefficients);
}
void
SteadyStateProblemBuilder::ConstructOperator()
{
  problem->eq_sys_operator =
      std::make_unique<hephaestus::EquationSystemOperator>(*(problem->pmesh),
                                                           problem->fespaces,
                                                           problem->gridfunctions,
                                                           problem->bc_map,
                                                           problem->coefficients,
                                                           problem->sources,
                                                           problem->solver_options);
  problem->eq_sys_operator->SetGridFunctions();
}

void
SteadyStateProblemBuilder::ConstructState()
{
  problem->F = new mfem::BlockVector(problem->eq_sys_operator->true_offsets); // Vector of dofs
  problem->eq_sys_operator->Init(*(problem->F)); // Set up initial conditions
}
} // namespace hephaestus

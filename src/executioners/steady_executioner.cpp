#include "steady_executioner.hpp"

namespace hephaestus {

SteadyExecutioner::SteadyExecutioner(const hephaestus::InputParameters &params)
    : Executioner(params),
      problem(params.GetParam<hephaestus::SteadyStateProblem *>("Problem")) {}

void SteadyExecutioner::Init() {
  // Set up DataCollections to track fields of interest.
  problem->outputs.SetGridFunctions(problem->gridfunctions);
  problem->outputs.Reset();
  problem->outputs.EnableGLVis(visualization);
}

void SteadyExecutioner::Solve() const {
  // Advance time step.
  problem->preprocessors.Solve();
  problem->GetOperator()->Solve(*(problem->F));
  problem->postprocessors.Solve();

  // Output data
  // Output timestep summary to console
  problem->outputs.Write();
}
void SteadyExecutioner::Execute() const { this->Solve(); }
} // namespace hephaestus

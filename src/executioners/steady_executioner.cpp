#include "steady_executioner.hpp"

namespace hephaestus {

SteadyExecutioner::SteadyExecutioner(const hephaestus::InputParameters &params)
    : Executioner(params),
      problem(
          params.GetParam<hephaestus::FrequencyDomainProblem *>("Problem")) {}

void SteadyExecutioner::Init() {
  problem->auxkernels.Solve();

  // Set up DataCollections to track fields of interest.
  for (auto const &[name, dc_] : problem->data_collections) {
    problem->outputs.RegisterOutputFields(dc_, problem->pmesh.get(),
                                          problem->gridfunctions);
    // Write initial fields to disk
    problem->outputs.WriteOutputFields(dc_, 0.0, 0);
  }

  // Initialize GLVis visualization and send the initial condition
  // by socket to a GLVis server.
  if (visualization) {
    problem->outputs.InitializeGLVis(problem->myid_, problem->gridfunctions);
    problem->outputs.DisplayToGLVis(problem->gridfunctions);
  }
}

void SteadyExecutioner::Solve() const {
  // Advance time step.
  problem->fd_operator.get()->Solve(*(problem->F));
  problem->auxkernels.Solve();

  // Output data
  // Output timestep summary to console
  problem->outputs.WriteConsoleSummary(problem->myid_, 1.0, 1);

  // Make sure all ranks have sent their 'v' solution before initiating
  // another set of GLVis connections (one from each rank):
  MPI_Barrier(problem->pmesh->GetComm());

  // Send output fields to GLVis for visualisation
  if (visualization) {
    problem->outputs.DisplayToGLVis(problem->gridfunctions);
  }

  // Save output fields at timestep to DataCollections
  for (auto const &[name, dc_] : problem->data_collections) {
    problem->outputs.WriteOutputFields(dc_, 1.0, 1);
  }
  problem->postprocessors.Update();
}
void SteadyExecutioner::Execute() const { this->Solve(); }
} // namespace hephaestus

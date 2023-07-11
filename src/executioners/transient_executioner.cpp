#include "transient_executioner.hpp"

namespace hephaestus {

TransientExecutioner::TransientExecutioner(
    const hephaestus::InputParameters &params)
    : Executioner(params),
      problem(params.GetParam<hephaestus::TimeDomainProblem *>("Problem")),
      t_step(params.GetParam<float>("TimeStep")),
      t_initial(params.GetParam<float>("StartTime")),
      t_final(params.GetParam<float>("EndTime")), t(t_initial), it(0),
      vis_steps(params.GetOptionalParam<int>("VisualisationSteps", 1)),
      last_step(false) {}

void TransientExecutioner::Init() {
  // Set up DataCollections to track fields of interest.
  for (auto const &[name, dc_] : problem->outputs.data_collections) {
    problem->outputs.RegisterOutputFields(dc_, problem->pmesh.get(),
                                          problem->gridfunctions);
    // Write initial fields to disk
    problem->outputs.WriteOutputFields(dc_, t, 0);
  }

  // Initialize GLVis visualization and send the initial condition
  // by socket to a GLVis server.
  if (visualization) {
    problem->outputs.InitializeGLVis(problem->myid_, problem->gridfunctions);
    problem->outputs.DisplayToGLVis(problem->gridfunctions);
  }
}

void TransientExecutioner::Step(double dt, int it) const {
  // Check if current time step is final
  if (t + dt >= t_final - dt / 2) {
    last_step = true;
  }

  // Advance time step.
  problem->preprocessors.Solve(t);
  problem->ode_solver->Step(*(problem->F), t, dt);
  problem->postprocessors.Solve(t);

  // Output data
  if (last_step || (it % vis_steps) == 0) {

    // Output timestep summary to console
    problem->outputs.WriteConsoleSummary(problem->myid_, t, it);

    // Make sure all ranks have sent their 'v' solution before initiating
    // another set of GLVis connections (one from each rank):
    MPI_Barrier(problem->pmesh->GetComm());

    // Send output fields to GLVis for visualisation
    if (visualization) {
      problem->outputs.DisplayToGLVis(problem->gridfunctions);
    }

    // Save output fields at timestep to DataCollections
    for (auto const &[name, dc_] : problem->outputs.data_collections) {
      problem->outputs.WriteOutputFields(dc_, t, it);
    }
  }
}

void TransientExecutioner::Solve() const {
  it++;
  Step(t_step, it);
}

void TransientExecutioner::Execute() const {
  // Initialise time gridfunctions
  t = t_initial;
  last_step = false;
  it = 0;
  while (last_step != true) {
    Solve();
  }
}

} // namespace hephaestus

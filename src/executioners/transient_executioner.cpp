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
  problem->outputs.SetGridFunctions(problem->gridfunctions);
  problem->outputs.Reset();
  problem->outputs.EnableGLVis(visualization);
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
    problem->outputs.Write(t);
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

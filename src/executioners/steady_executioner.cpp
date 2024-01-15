#include "steady_executioner.hpp"

namespace hephaestus
{

SteadyExecutioner::SteadyExecutioner(const hephaestus::InputParameters & params)
  : Executioner(params), problem(params.GetParam<hephaestus::SteadyStateProblem *>("Problem"))
{
}

void
SteadyExecutioner::Solve() const
{
  // Advance time step.
  problem->preprocessors.Solve();
  problem->GetOperator()->Solve(*(problem->F));
  problem->postprocessors.Solve();

  // Output data
  // Output timestep summary to console
  problem->outputs.Write();
}
void
SteadyExecutioner::Execute() const
{
  Solve();
}
} // namespace hephaestus

#include "auxsolvers.hpp"

namespace hephaestus
{

void
AuxSolvers::Init(const hephaestus::GridFunctions & gridfunctions,
                 hephaestus::Coefficients & coefficients)
{
  for (const auto & [name, auxsolver] : *this)
  {
    logger.info("Initialising {} AuxSolver", name);
    spdlog::stopwatch sw;
    auxsolver->Init(gridfunctions, coefficients);
    _aux_queue.push_back(auxsolver);
    logger.info("{} Init: {} seconds", name, sw);
  }

  std::sort(_aux_queue.begin(), _aux_queue.end(), AuxSolver::PriorityComparator);
}

void
AuxSolvers::Solve(double t)
{
  for (auto & auxsolver : _aux_queue)
  {
    spdlog::stopwatch sw;
    auxsolver->Solve(t);
    logger.info("AuxSolver Solve: {} seconds", sw);
  }
}

} // namespace hephaestus

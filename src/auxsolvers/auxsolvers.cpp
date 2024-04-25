#include "auxsolvers.hpp"

namespace hephaestus
{

void
AuxSolvers::Init(const hephaestus::GridFunctions & gridfunctions,
                 hephaestus::Coefficients & coefficients)
{
  _aux_queue.clear();

  for (const auto & [name, auxsolver] : *this)
  {
    logger.info("Initialising {} AuxSolver", name);
    spdlog::stopwatch sw;
    auxsolver->Init(gridfunctions, coefficients);
    _aux_queue.emplace_back(std::pair(auxsolver, name));
    logger.info("{} Init: {} seconds", name, sw);
  }

  std::sort(_aux_queue.begin(), _aux_queue.end(), AuxSolver::PriorityComparator);
}

void
AuxSolvers::Update(const hephaestus::GridFunctions & gridfunctions,
                   hephaestus::Coefficients & coefficients)
{
  _aux_queue.clear();

  for ([[maybe_unused]] const auto & [name, auxsolver] : *this)
  {
    logger.debug("Update called for auxsolver '{}'.", name);

    auxsolver->Init(gridfunctions, coefficients);

    auto pair = std::make_pair(auxsolver, name);
    _aux_queue.emplace_back(auxsolver, name);
  }
}

void
AuxSolvers::Solve(double t)
{
  for (auto & aux_pair : _aux_queue)
  {
    spdlog::stopwatch sw;
    aux_pair.first->Solve(t);
    logger.info("{} AuxSolver Solve: {} seconds", aux_pair.second, sw);
  }
}

} // namespace hephaestus

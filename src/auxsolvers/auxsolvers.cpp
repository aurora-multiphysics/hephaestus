#include "auxsolvers.hpp"

namespace hephaestus
{

void
AuxSolvers::Init(const hephaestus::GridFunctions & gridfunctions,
                 hephaestus::Coefficients & coefficients)
{

  for (const auto & [name, auxsolver] : *this)
  {
    auxsolver->Init(gridfunctions, coefficients);
    _aux_queue.push_back(auxsolver);
  }

  std::sort(_aux_queue.begin(), _aux_queue.end(), AuxSolver::PriorityComparator);
}

void
AuxSolvers::Solve(double t)
{
  for (auto & auxsolver : _aux_queue)
  {
    auxsolver->Solve(t);
  }
}

} // namespace hephaestus

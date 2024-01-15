#include "auxsolvers.hpp"

namespace hephaestus
{

void
AuxSolvers::Init(const hephaestus::GridFunctions & gridfunctions,
                 hephaestus::Coefficients & coefficients)
{

  for (const auto & [name, auxsolver] : GetMap())
  {
    auxsolver->Init(gridfunctions, coefficients);
    aux_queue.push_back(auxsolver);
  }

  std::sort(aux_queue.begin(), aux_queue.end(), AuxSolver::comparator);
}

void
AuxSolvers::Solve(double t)
{
  for (auto & auxsolver : aux_queue)
  {
    auxsolver->Solve(t);
  }
}

} // namespace hephaestus

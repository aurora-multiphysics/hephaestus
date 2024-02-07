#include "mesh_delegate.hpp"

namespace hephaestus
{
void
MeshDelegator::MeshDidUpdate()
{
  for (const auto & delegate : _delegates)
  {
    delegate->OnMeshDidUpdate();
  }
}
}

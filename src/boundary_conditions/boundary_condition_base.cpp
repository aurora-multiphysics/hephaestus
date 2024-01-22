#include "boundary_condition_base.hpp"

#include <utility>

namespace hephaestus
{

BoundaryCondition::BoundaryCondition(std::string name_, mfem::Array<int> bdr_attributes_)
  : name(std::move(name_)), bdr_attributes(std::move(bdr_attributes_))
{
}

mfem::Array<int>
BoundaryCondition::GetMarkers(mfem::Mesh & mesh)
{
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes, markers);
  return markers;
}

} // namespace hephaestus

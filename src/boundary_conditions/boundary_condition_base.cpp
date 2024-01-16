#include "boundary_condition_base.hpp"

namespace hephaestus
{

BoundaryCondition::BoundaryCondition(const std::string & name_, mfem::Array<int> bdr_attributes_)
  : name(name_), bdr_attributes(bdr_attributes_)
{
}

mfem::Array<int>
BoundaryCondition::getMarkers(mfem::Mesh & mesh)
{
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes, markers);
  return markers;
}

} // namespace hephaestus

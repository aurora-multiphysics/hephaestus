#include "boundary_condition_base.hpp"

#include <utility>

namespace hephaestus
{

BoundaryCondition::BoundaryCondition(std::string name_, mfem::Array<int> bdr_attributes_)
  : _name(std::move(name_)), _bdr_attributes(std::move(bdr_attributes_))
{
}

mfem::Array<int>
BoundaryCondition::GetMarkers(mfem::Mesh & mesh)
{
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), _bdr_attributes, _markers);
  return _markers;
}

void
BoundaryCondition::Update(mfem::Mesh & mesh)
{
  // If the mesh has been refined, the number of boundary attributes will have
  // changed! We will need to map from our old bdr_attributes array to the new
  // one.
  if (mesh.bdr_attributes.Size() == _bdr_attributes.Size())
  {
    return;
  }

  // Ensure that we have non-empty arrays.
  if (mesh.bdr_attributes.Size() == 0 || _bdr_attributes.Size() == 0)
  {
    MFEM_ABORT("Encountered empty boundary attribute arrays.");
  }

  // NB: assume x2 for 2D uniform refinement; x4 for 3D.
  int refinement_factor = (mesh.bdr_attributes.Size() / _bdr_attributes.Size());

  if (refinement_factor != 2 && refinement_factor != 4)
  {
    MFEM_ABORT("Cannot handle refinement factor: " << refinement_factor);
  }

  mfem::Array<int> new_bdr_attributes(mesh.bdr_attributes.Size());

  for (int i = 0; i < _bdr_attributes.Size(); i++)
  {
    int attr = _bdr_attributes[i];

    for (int k = 0; k < refinement_factor; k++)
    {
      new_bdr_attributes[i * refinement_factor + k] = attr;
    }
  }

  _bdr_attributes = new_bdr_attributes;
}

} // namespace hephaestus

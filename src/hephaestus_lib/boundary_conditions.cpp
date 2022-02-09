#include "boundary_conditions.hpp"


BCMap::BCMap(const std::string & bc_name,
    Array<int> bdr_attr)
    :name(bc_name),
    bdr_attributes(bdr_attr)
{
}

Array<int> BCMap::getMarkers(Mesh & mesh)
{
    mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes, markers);
    return markers;
}
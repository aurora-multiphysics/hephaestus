#include "boundary_conditions.hpp"

namespace hephaestus
{

BoundaryCondition::BoundaryCondition()
{
}

BoundaryCondition::BoundaryCondition(const std::string & boundary_name,
    mfem::Array<int> boundary_ids)
    :name(boundary_name),
    bdr_attributes(boundary_ids)
{
}

mfem::Array<int> BoundaryCondition::getMarkers(mfem::Mesh & mesh)
{
    std::cout<<"MAX ATTRIBUTES"<<mesh.bdr_attributes.Max();
    mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes, markers);
    return markers;
}

BCMap::BCMap()
{
}

void BCMap::setBC(std::string bc_name, BoundaryCondition bc)
{
    bc_map.insert(std::pair<std::string, BoundaryCondition *>(bc_name, new BoundaryCondition(bc)));
}

BoundaryCondition BCMap::getBC(std::string bc_name)
{
    BoundaryCondition * test_bc = bc_map[bc_name];
    return *test_bc;
}

NeumannBC::NeumannBC()
{
}

NeumannBC::NeumannBC(const std::string & boundary_name, mfem::Array<int> boundary_ids)
:BoundaryCondition(boundary_name, boundary_ids)
{
}

} //hephaestus
#include "inputs.hpp"

namespace hephaestus{

Inputs::Inputs(const std::string & mesh_file,
        const int order,
        const BCMap & boundary_conditions,
        const MaterialMap & materials)
    :
    _mesh_file(mesh_file),
    _order(order),
    bc_map(boundary_conditions),
    material_map(materials)
{
}

} //namespace hephaestus
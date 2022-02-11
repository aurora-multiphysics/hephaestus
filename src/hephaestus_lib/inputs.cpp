#include "inputs.hpp"

namespace hephaestus{

Inputs::Inputs(const std::string & problem_type,
        const std::string & mesh_file,
        const BCMap & boundary_conditions,
        const MaterialMap & materials)
    :
    _problem_type(problem_type),
    _mesh_file(mesh_file),
    bc_map(boundary_conditions),
    material_map(materials)
{
}

} //namespace hephaestus
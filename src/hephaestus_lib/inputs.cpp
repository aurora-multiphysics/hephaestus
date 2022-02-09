#include "inputs.hpp"


Inputs::Inputs(const std::string & problem_type,
        const std::string & mesh_file,
        const std::vector<BCMap> & boundary_conditions)
    :
    _problem_type(problem_type),
    _mesh_file(mesh_file),
    bc_maps(boundary_conditions)
{
}


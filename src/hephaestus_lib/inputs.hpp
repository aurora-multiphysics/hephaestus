#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "boundary_conditions.hpp"
#include "materials.hpp"
#include "joule_solver.hpp"

namespace hephaestus{

class Inputs
{
    public:
    Inputs(const std::string & mesh_file,
        const int order,
        const BCMap & boundary_conditions,
        const MaterialMap & materials);

    std::string _mesh_file;
    int _order;
    BCMap bc_map;
    MaterialMap material_map;
};

} //namespace hephaestus
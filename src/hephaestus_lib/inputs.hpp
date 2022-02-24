#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "boundary_conditions.hpp"
#include "materials.hpp"
#include "executioner.hpp"
#include "joule_solver.hpp"

namespace hephaestus{

class Inputs
{
    public:
    Inputs(){};
    Inputs(const std::string & mesh_file,
        const std::string & formulation,
        const int order,
        const BCMap & boundary_conditions,
        const MaterialMap & materials,
        const Executioner & executioner_);

    std::string _mesh_file;
    std::string _formulation;
    int _order;
    BCMap bc_map;
    MaterialMap material_map;
    Executioner executioner;
};

} //namespace hephaestus
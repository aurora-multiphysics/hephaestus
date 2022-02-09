#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "boundary_conditions.hpp"
#include "joule_solver.hpp"

class Inputs
{
    public:
    Inputs(const std::string & problem_type,
        const std::string & mesh_file,
        const std::vector<BCMap> & boundary_conditions);

    std::string _problem_type;
    std::string _mesh_file;
    std::vector<BCMap> bc_maps;
};
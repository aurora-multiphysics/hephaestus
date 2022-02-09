#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "boundary_conditions.hpp"
#include "joule_solver.hpp"

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;

class Inputs
{
    public:
    Inputs(const std::vector<BCMap> & boundary_conditions);

    std::vector<BCMap> bc_maps;
};
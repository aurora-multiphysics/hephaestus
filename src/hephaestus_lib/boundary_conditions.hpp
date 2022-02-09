#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "joule_solver.hpp"
#include "mesh_extras.hpp"

class BCMap
{
    public:
    BCMap(const std::string & bc_name,
          mfem::Array<int> boundary_ids);
    mfem::Array<int> getMarkers(mfem::Mesh & mesh);

    std::string name;
    mfem::Array<int> bdr_attributes;
    mfem::Array<int> markers;
};
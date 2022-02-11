#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "joule_solver.hpp"
#include "mesh_extras.hpp"


namespace hephaestus
{

class BoundaryCondition
{
    public:
    BoundaryCondition();
    BoundaryCondition(const std::string & boundary_name,
          mfem::Array<int> boundary_ids);
    mfem::Array<int> getMarkers(mfem::Mesh & mesh);

    std::string name;
    mfem::Array<int> bdr_attributes;
    mfem::Array<int> markers;
};

class BCMap
{
    public:
    BCMap();

    void setBC(std::string bc_name, BoundaryCondition bc);
    BoundaryCondition getBC(std::string bc_name);

    std::map<std::string, hephaestus::BoundaryCondition> bc_map;
};

}
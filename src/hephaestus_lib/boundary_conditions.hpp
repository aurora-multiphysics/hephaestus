#pragma once
#include <memory>
#include <iostream>
#include <fstream>
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
    std::function<double(const mfem::Vector&, double)> scalar_func;
    std::function<void(const mfem::Vector&, double, mfem::Vector&)> vector_func;
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
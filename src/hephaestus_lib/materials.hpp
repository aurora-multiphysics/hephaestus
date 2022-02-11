#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "joule_solver.hpp"
#include "mesh_extras.hpp"

namespace hephaestus
{

class Material
{
    public:
    Material(const std::string & material_name,
             int material_block_id);

    std::string name;
    int block_id;
    std::map<std::string, double> properties;


    void setMaterialProperty(std::string property_name, double property_value);
    double getMaterialProperty(std::string property_name);
};

class MaterialMap
{
    public:
    MaterialMap(std::vector<Material> mats);

    std::map<int, double> getBlockPropertyMap(std::string property_name);

    std::vector<Material> materials;
};

} // namespace hephaestus
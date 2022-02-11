#include "materials.hpp"

namespace hephaestus
{

Material::Material(const std::string & material_name,
        int material_block_id)
        :name(material_name),
        block_id(material_block_id)
{
}

void Material::setMaterialProperty(std::string property_name, double property_value)
{
    properties.insert(std::pair<std::string, double>(property_name, property_value));
}

double Material::getMaterialProperty(std::string property_name)
{
    return properties[property_name];
}

MaterialMap::MaterialMap()
{
}

MaterialMap::MaterialMap(std::vector<Material> mats)
    :materials(mats)
{
}

std::map<int, double> MaterialMap::getBlockPropertyMap(std::string property_name)
{
    std::map<int, double> block_property_map;
    int block_id;
    double property_value;
    for (std::size_t i = 0; i < materials.size(); i++){
        Material material(materials[i]);
        block_id = material.block_id;
        property_value = material.properties[property_name];
        block_property_map.insert(std::pair<int, double>(block_id,property_value));
   }
    return block_property_map;
}

} // namespace hephaestus
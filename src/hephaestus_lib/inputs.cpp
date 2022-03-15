#include "inputs.hpp"

namespace hephaestus {
Inputs::Inputs(const std::string& mesh_file, const std::string& formulation,
               const int order, const BCMap& boundary_conditions,
               const MaterialMap& materials, const Executioner& executioner_)
    : _mesh_file(mesh_file),
      _formulation(formulation),
      _order(order),
      bc_map(boundary_conditions),
      material_map(materials),
      executioner(executioner_) {}

}  // namespace hephaestus
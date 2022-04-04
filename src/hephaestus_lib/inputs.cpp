#include "inputs.hpp"

namespace hephaestus {
Inputs::Inputs(const std::string &mesh_file, const std::string &formulation,
               const int order, const BCMap &boundary_conditions,
               const DomainProperties &domain_properties_,
               const Executioner &executioner_)
    : _mesh_file(mesh_file), _formulation(formulation), _order(order),
      bc_map(boundary_conditions), domain_properties(domain_properties_),
      executioner(executioner_) {}

} // namespace hephaestus
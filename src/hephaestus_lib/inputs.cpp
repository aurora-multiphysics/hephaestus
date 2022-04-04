#include "inputs.hpp"

namespace hephaestus {
Inputs::Inputs(const std::string &mesh_file_, const std::string &formulation_,
               const int order_, const BCMap &bc_map_,
               const DomainProperties &domain_properties_,
               const Executioner &executioner_)
    : mesh_file(mesh_file_), formulation(formulation_), order(order_),
      bc_map(bc_map_), domain_properties(domain_properties_),
      executioner(executioner_) {}

} // namespace hephaestus
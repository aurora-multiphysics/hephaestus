#include "inputs.hpp"

namespace hephaestus {
Inputs::Inputs(const mfem::Mesh &mesh_, const std::string &formulation_,
               const int order_, const BCMap &bc_map_,
               const DomainProperties &domain_properties_,
               const Executioner &executioner_, Outputs outputs_)
    : mesh(mesh_), formulation(formulation_), order(order_), bc_map(bc_map_),
      domain_properties(domain_properties_), executioner(executioner_),
      outputs(outputs_) {}

} // namespace hephaestus
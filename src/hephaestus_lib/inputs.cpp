#include "inputs.hpp"

namespace hephaestus {
Executioner::Executioner(const std::string &type_, const double dt_,
                         const double t_initial_, const double t_final_)
    : type(type_), dt(dt_), t_initial(t_initial_), t_final(t_final_) {}

Inputs::Inputs(const mfem::Mesh &mesh_, const std::string &formulation_,
               const int order_, const BCMap &bc_map_,
               const DomainProperties &domain_properties_,
               const Executioner &executioner_, Outputs outputs_)
    : mesh(mesh_), formulation(formulation_), order(order_), bc_map(bc_map_),
      domain_properties(domain_properties_), executioner(executioner_),
      outputs(outputs_) {}

} // namespace hephaestus

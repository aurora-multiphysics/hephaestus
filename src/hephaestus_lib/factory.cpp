#include "factory.hpp"

namespace hephaestus {

hephaestus::TransientFormulation *
FormulationFactory::createTransientFormulation(
    std::string &formulation, mfem::ParMesh &pmesh, int order,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  if (formulation == "EBForm") {
    return new hephaestus::EBDualSolver(pmesh, order, bc_map,
                                        domain_properties);
  } else if (formulation == "HJForm") {
    return new hephaestus::HJDualSolver(pmesh, order, bc_map,
                                        domain_properties);
  } else if (formulation == "HForm") {
    return new hephaestus::HFormSolver(pmesh, order, bc_map, domain_properties);
  }
  return nullptr;
}

} // namespace hephaestus

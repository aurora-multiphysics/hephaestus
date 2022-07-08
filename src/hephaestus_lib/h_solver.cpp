#include "h_solver.hpp"

namespace hephaestus {

HSolver::HSolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
                 hephaestus::DomainProperties &domain_properties)
    : HCurlSolver(pmesh, order, bc_map, domain_properties) {}

void HSolver::SetVariableNames() {
  p_name = "magnetic_potential";
  p_display_name = "Scalar Potential (V)";

  u_name = "magnetic_field";
  u_display_name = "Magnetic Field (H)";

  v_name = "current_density";
  v_display_name = "Current Density (J)";
}

void HSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef,
      domain_properties.scalar_property_map["electrical_conductivity"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["magnetic_permeability"];
}

} // namespace hephaestus

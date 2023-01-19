#include "eb_dual_solver.hpp"

namespace hephaestus {

EBDualSolver::EBDualSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : DualSolver(pmesh, order, fespaces, variables, bc_map, domain_properties,
                 sources, solver_options) {}

void EBDualSolver::RegisterVariables() {
  u_name = "electric_field";
  u_display_name = "Electric Field (E)";

  v_name = "magnetic_flux_density";
  v_display_name = "Magnetic Flux Density (B)";

  _variables.Register(u_name, &u_, false);
  _variables.Register(v_name, &v_, false);
}

void EBDualSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
      0) {
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("magnetic_permeability")));
  }
  if (domain_properties.scalar_property_map.count("electrical_conductivity") ==
      0) {
    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("electrical_conductivity")));
  }

  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef, domain_properties.scalar_property_map["magnetic_permeability"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
}

} // namespace hephaestus

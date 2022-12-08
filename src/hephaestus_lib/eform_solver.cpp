#include "eform_solver.hpp"

namespace hephaestus {

EFormSolver::EFormSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources)
    : HCurlSolver(pmesh, order, fespaces, variables, bc_map, domain_properties,
                  sources) {}

void EFormSolver::SetVariableNames() {
  p_name = "electric_potential";
  p_display_name = "Electric Potential (V)";

  u_name = "electric_field";
  u_display_name = "Electric Field (E)";
}

void EFormSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("electrical_conductivity") ==
      0) {
    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("electrical_conductivity")));
  }
  if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
      0) {
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("magnetic_permeability")));
  }
  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef, domain_properties.scalar_property_map["magnetic_permeability"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
}

} // namespace hephaestus

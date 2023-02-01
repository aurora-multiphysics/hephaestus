#include "aform_solver.hpp"

namespace hephaestus {

AFormSolver::AFormSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : HCurlSolver(pmesh, order, fespaces, variables, bc_map, domain_properties,
                  sources, solver_options) {}

void AFormSolver::RegisterVariables() {
  u_name = "magnetic_vector_potential";
  u_display_name = "Magnetic Vector Potential (A)";
  curl_u_name = "magnetic_flux_density";

  HCurlSolver::RegisterVariables();
}

void AFormSolver::SetMaterialCoefficients(
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

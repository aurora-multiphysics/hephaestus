#include "hform_solver.hpp"

namespace hephaestus {

HFormSolver::HFormSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Kernels &kernels, hephaestus::InputParameters &solver_options)
    : HCurlSolver(pmesh, order, fespaces, variables, bc_map, domain_properties,
                  kernels, solver_options) {}

void HFormSolver::RegisterVariables() {
  u_name = "magnetic_field";
  u_display_name = "Magnetic Field (H)";

  _variables.Register(u_name, &u_, false);
  _variables.Register("current_density", &curl_u_, false);
}

void HFormSolver::SetMaterialCoefficients(
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

  alpha_coef_name = std::string("electrical_resistivity");
  beta_coef_name = std::string("magnetic_permeability");

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["electrical_conductivity"],
          fracFunc);

  alphaCoef = domain_properties.scalar_property_map[alpha_coef_name];
  betaCoef = domain_properties.scalar_property_map[beta_coef_name];
}

} // namespace hephaestus

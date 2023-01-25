#include "eform_solver.hpp"

namespace hephaestus {

EFormSolver::EFormSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : HCurlSolver(pmesh, order, fespaces, variables, bc_map, domain_properties,
                  sources, solver_options) {}

void EFormSolver::RegisterVariables() {
  u_name = "electric_field";
  u_display_name = "Electric Field (E)";

  _variables.Register(u_name, &u_, false);
  _variables.Register("-dB_dt", &curl_u_, false);
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
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);

  alphaCoef = domain_properties.scalar_property_map[alpha_coef_name];
  betaCoef = domain_properties.scalar_property_map[beta_coef_name];
}

} // namespace hephaestus

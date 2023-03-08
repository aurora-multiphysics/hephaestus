#include "eb_dual_solver.hpp"

namespace hephaestus {

EBDualFormulation::EBDualFormulation() : DualFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  CreateEquationSystem();
}

hephaestus::TimeDependentEquationSystem *
EBDualFormulation::CreateEquationSystem() {
  std::vector<std::string> state_var_names;
  state_var_names.resize(2);
  state_var_names.at(0) = "electric_field";
  state_var_names.at(1) = "magnetic_flux_density";

  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("VariableNames", state_var_names);
  weak_form_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  equation_system = new hephaestus::CurlCurlEquationSystem(weak_form_params);
  return equation_system;
}

hephaestus::TimeDomainEquationSystemOperator *
EBDualFormulation::CreateTimeDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  td_operator = new hephaestus::EBDualOperator(
      pmesh, order, fespaces, variables, bc_map, domain_properties, sources,
      solver_options);
  td_operator->SetEquationSystem(equation_system);
  return td_operator;
};

void EBDualFormulation::RegisterCoefficients(
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

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);
}

EBDualOperator::EBDualOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : DualOperator(pmesh, order, fespaces, variables, bc_map, domain_properties,
                   sources, solver_options) {}

void EBDualOperator::RegisterVariables() {
  u_name = "electric_field";
  u_display_name = "Electric Field (E)";

  v_name = "magnetic_flux_density";
  v_display_name = "Magnetic Flux Density (B)";

  _variables.Register(u_name, &u_, false);
  _variables.Register(v_name, &v_, false);
}

void EBDualOperator::SetMaterialCoefficients(
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

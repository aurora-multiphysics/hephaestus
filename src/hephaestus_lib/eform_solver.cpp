#include "eform_solver.hpp"

namespace hephaestus {

EFormulation::EFormulation() : HCurlFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  CreateEquationSystem();
}

hephaestus::TimeDependentEquationSystem *EFormulation::CreateEquationSystem() {
  std::vector<std::string> state_var_names;
  state_var_names.resize(1);
  state_var_names.at(0) = "electric_field";

  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("VariableNames", state_var_names);
  weak_form_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  equation_system = new hephaestus::CurlCurlEquationSystem(weak_form_params);
  return equation_system;
}

void EFormulation::RegisterCoefficients(
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

// EFormSolver::EFormSolver(
//     mfem::ParMesh &pmesh, int order,
//     mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//     mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//     hephaestus::BCMap &bc_map, hephaestus::DomainProperties
//     &domain_properties, hephaestus::Sources &sources,
//     hephaestus::InputParameters &solver_options) : HCurlSolver(pmesh, order,
//     fespaces, variables, bc_map, domain_properties,
//                   sources, solver_options) {

//   state_var_names.resize(1);
//   state_var_names.at(0) = "electric_field";

//   aux_var_names.resize(1);
//   aux_var_names.at(0) = "-dB_dt";
// }

// void EFormSolver::SetMaterialCoefficients(
//     hephaestus::DomainProperties &domain_properties) {
//   if (domain_properties.scalar_property_map.count("electrical_conductivity")
//   ==
//       0) {
//     domain_properties.scalar_property_map["electrical_conductivity"] =
//         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
//             std::string("electrical_conductivity")));
//   }
//   if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
//       0) {
//     domain_properties.scalar_property_map["magnetic_permeability"] =
//         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
//             std::string("magnetic_permeability")));
//   }

//   alpha_coef_name = std::string("magnetic_reluctivity");
//   beta_coef_name = std::string("electrical_conductivity");

//   domain_properties.scalar_property_map[alpha_coef_name] =
//       new mfem::TransformedCoefficient(
//           &oneCoef,
//           domain_properties.scalar_property_map["magnetic_permeability"],
//           fracFunc);
// }

} // namespace hephaestus

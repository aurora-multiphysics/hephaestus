#include "hj_dual_solver.hpp"

namespace hephaestus {

HJDualFormulation::HJDualFormulation() : DualFormulation() {
  alpha_coef_name = std::string("electrical_resistivity");
  beta_coef_name = std::string("magnetic_permeability");
  CreateEquationSystem();
}

hephaestus::TimeDependentEquationSystem *
HJDualFormulation::CreateEquationSystem() {
  std::vector<std::string> state_var_names;
  state_var_names.resize(2);
  state_var_names.at(0) = "magnetic_field";
  state_var_names.at(1) = "current_density";

  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("VariableNames", state_var_names);
  weak_form_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  equation_system = new hephaestus::CurlCurlEquationSystem(weak_form_params);
  return equation_system;
}

void HJDualFormulation::RegisterCoefficients(
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
          domain_properties.scalar_property_map["electrical_conductivity"],
          fracFunc);
}

// HJDualSolver::HJDualSolver(
//     mfem::ParMesh &pmesh, int order,
//     mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//     mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//     hephaestus::BCMap &bc_map, hephaestus::DomainProperties
//     &domain_properties, hephaestus::Sources &sources,
//     hephaestus::InputParameters &solver_options) : DualSolver(pmesh, order,
//     fespaces, variables, bc_map, domain_properties,
//                  sources, solver_options) {}

// void HJDualSolver::RegisterVariables() {
//   u_name = "magnetic_field";
//   u_display_name = "Magnetic Field (H)";

//   v_name = "current_density";
//   v_display_name = "Current Density (J)";

//   _variables.Register(u_name, &u_, false);
//   _variables.Register(v_name, &v_, false);
// }

// void HJDualSolver::SetMaterialCoefficients(
//     hephaestus::DomainProperties &domain_properties) {
//   if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
//       0) {
//     domain_properties.scalar_property_map["magnetic_permeability"] =
//         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
//             std::string("magnetic_permeability")));
//   }
//   if (domain_properties.scalar_property_map.count("electrical_conductivity")
//   ==
//       0) {
//     domain_properties.scalar_property_map["electrical_conductivity"] =
//         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
//             std::string("electrical_conductivity")));
//   }
//   alphaCoef = new mfem::TransformedCoefficient(
//       &oneCoef,
//       domain_properties.scalar_property_map["electrical_conductivity"],
//       fracFunc);
//   betaCoef = domain_properties.scalar_property_map["magnetic_permeability"];
// }

} // namespace hephaestus

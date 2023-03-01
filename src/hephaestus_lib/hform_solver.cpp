#include "hform_solver.hpp"

namespace hephaestus {

HFormulation::HFormulation() : HCurlFormulation() {
  alpha_coef_name = std::string("electrical_resistivity");
  beta_coef_name = std::string("magnetic_permeability");
  CreateEquationSystem();
}

hephaestus::TimeDependentEquationSystem *HFormulation::CreateEquationSystem() {
  std::vector<std::string> state_var_names;
  state_var_names.resize(1);
  state_var_names.at(0) = "magnetic_field";

  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("VariableNames", state_var_names);
  weak_form_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  equation_system = new hephaestus::CurlCurlEquationSystem(weak_form_params);
  return equation_system;
}

void HFormulation::RegisterAuxKernels(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::AuxKernels &auxkernels) {
  std::vector<std::string> aux_var_names;
  std::string j_field_name = "current_density";
  if (variables.Get(j_field_name) != NULL) {
    // if (myid_ == 0) {
    //   std::cout << j_field_name << " found in variables: building auxvar "
    //             << std::endl;
    // }
    hephaestus::InputParameters j_field_aux_params;
    j_field_aux_params.SetParam("VariableName", std::string("magnetic_field"));
    j_field_aux_params.SetParam("CurlVariableName", j_field_name);
    auxkernels.Register("_magnetic_flux_density_aux",
                        new hephaestus::CurlAuxKernel(j_field_aux_params),
                        true);
  }
}

void HFormulation::RegisterCoefficients(
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

// HFormSolver::HFormSolver(
//     mfem::ParMesh &pmesh, int order,
//     mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//     mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//     hephaestus::BCMap &bc_map, hephaestus::DomainProperties
//     &domain_properties, hephaestus::Sources &sources,
//     hephaestus::InputParameters &solver_options) : HCurlSolver(pmesh, order,
//     fespaces, variables, bc_map, domain_properties,
//                   sources, solver_options) {

//   state_var_names.resize(1);
//   state_var_names.at(0) = "magnetic_field";

//   aux_var_names.resize(1);
//   aux_var_names.at(0) = "current_density";
// }

// void HFormSolver::RegisterAuxKernels(hephaestus::AuxKernels &auxkernels) {
//   for (auto &active_aux_var_name : active_aux_var_names) {
//     // Check if current density should be added as auxvar
//     if (active_aux_var_name == aux_var_names.at(0)) {
//       if (myid_ == 0) {
//         std::cout << active_aux_var_name
//                   << " found in variables: building auxvar " << std::endl;
//       }
//       hephaestus::InputParameters current_density_aux_params;
//       current_density_aux_params.SetParam("VariableName",
//                                           state_var_names.at(0));
//       current_density_aux_params.SetParam("CurlVariableName",
//                                           active_aux_var_name);
//       auxkernels.Register(
//           "_current_density_aux",
//           new hephaestus::CurlAuxKernel(current_density_aux_params), true);
//     }
//   }
// }
// void HFormSolver::SetMaterialCoefficients(
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
//   alpha_coef_name = std::string("electrical_resistivity");
//   beta_coef_name = std::string("magnetic_permeability");

//   domain_properties.scalar_property_map[alpha_coef_name] =
//       new mfem::TransformedCoefficient(
//           &oneCoef,
//           domain_properties.scalar_property_map["electrical_conductivity"],
//           fracFunc);
// }

} // namespace hephaestus

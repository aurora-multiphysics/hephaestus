#include "hform_solver.hpp"

namespace hephaestus {

HFormSolver::HFormSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : HCurlSolver(pmesh, order, fespaces, variables, bc_map, domain_properties,
                  sources, solver_options) {

  state_var_names.resize(1);
  state_var_names[0] = "magnetic_field";

  aux_var_names.resize(1);
  aux_var_names[0] = "current_density";
}

void HFormSolver::RegisterAuxKernels(hephaestus::AuxKernels &auxkernels) {
  hephaestus::InputParameters current_density_aux_params;
  current_density_aux_params.SetParam("VariableName", u_name);
  current_density_aux_params.SetParam("CurlVariableName", curl_u_name);
  auxkernels.Register("_current_density_aux",
                      new hephaestus::CurlAuxKernel(current_density_aux_params),
                      true);
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
  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef,
      domain_properties.scalar_property_map["electrical_conductivity"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["magnetic_permeability"];
}

} // namespace hephaestus

#include "hj_dual_formulation.hpp"

namespace hephaestus {

HJDualFormulation::HJDualFormulation() : DualFormulation() {
  alpha_coef_name = std::string("electrical_resistivity");
  beta_coef_name = std::string("magnetic_permeability");
  h_curl_var_name = std::string("magnetic_field");
  h_div_var_name = std::string("current_density");
}

void HJDualFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &domain_properties =
      this->GetProblem()->domain_properties;
  if (!domain_properties.scalar_property_map.Has("magnetic_permeability")) {
    MFEM_ABORT("magnetic_permeability coefficient not found.");
  }
  if (!domain_properties.scalar_property_map.Has("electrical_conductivity")) {
    MFEM_ABORT("electrical_conductivity coefficient not found.");
  }
  domain_properties.scalar_property_map.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map.Get("electrical_conductivity"),
          fracFunc),
      true);
}

} // namespace hephaestus

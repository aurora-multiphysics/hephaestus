#include "eb_dual_formulation.hpp"

namespace hephaestus {

EBDualFormulation::EBDualFormulation() : DualFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  h_curl_var_name = std::string("electric_field");
  h_div_var_name = std::string("magnetic_flux_density");
}

void EBDualFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &domain_properties =
      this->GetProblem()->domain_properties;

  if (!domain_properties.scalar_property_map.Has("magnetic_permeability")) {
    MFEM_ABORT("magnetic_permeability coefficient not found.");
  }
  if (!domain_properties.scalar_property_map.Has(beta_coef_name)) {
    MFEM_ABORT(beta_coef_name + " coefficient not found.");
  }

  domain_properties.scalar_property_map.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map.Get("magnetic_permeability"),
          fracFunc),
      true);
}

} // namespace hephaestus

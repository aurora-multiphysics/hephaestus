#include "eb_dual_formulation.hpp"

namespace hephaestus {

EBDualFormulation::EBDualFormulation() : DualFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  h_curl_var_name = std::string("electric_field");
  h_div_var_name = std::string("magnetic_flux_density");
}

void EBDualFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;

  if (!coefficients.scalars.Has("magnetic_permeability")) {
    MFEM_ABORT("magnetic_permeability coefficient not found.");
  }
  if (!coefficients.scalars.Has(beta_coef_name)) {
    MFEM_ABORT(beta_coef_name + " coefficient not found.");
  }

  coefficients.scalars.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get("magnetic_permeability"),
          fracFunc),
      true);
}

} // namespace hephaestus

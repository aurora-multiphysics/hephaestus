#include "hj_dual_formulation.hpp"

namespace hephaestus {

HJDualFormulation::HJDualFormulation() : DualFormulation() {
  alpha_coef_name = std::string("electrical_resistivity");
  beta_coef_name = std::string("magnetic_permeability");
  h_curl_var_name = std::string("magnetic_field");
  h_div_var_name = std::string("current_density");
}

void HJDualFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has("magnetic_permeability")) {
    MFEM_ABORT("magnetic_permeability coefficient not found.");
  }
  if (!coefficients.scalars.Has("electrical_conductivity")) {
    MFEM_ABORT("electrical_conductivity coefficient not found.");
  }
  coefficients.scalars.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get("electrical_conductivity"),
          fracFunc),
      true);
}

} // namespace hephaestus

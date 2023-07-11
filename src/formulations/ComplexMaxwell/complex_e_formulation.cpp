#include "complex_e_formulation.hpp"

namespace hephaestus {

ComplexEFormulation::ComplexEFormulation()
    : hephaestus::ComplexMaxwellFormulation() {
  h_curl_var_name = std::string("electric_field");

  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  zeta_coef_name = std::string("dielectric_permittivity");
};

} // namespace hephaestus

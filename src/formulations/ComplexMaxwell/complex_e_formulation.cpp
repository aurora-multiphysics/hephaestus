#include "complex_e_formulation.hpp"

namespace hephaestus {

ComplexEFormulation::ComplexEFormulation(
    const std::string &magnetic_reluctivity_name,
    const std::string &electric_conductivity_name,
    const std::string &dielectric_permittivity_name,
    const std::string &frequency_coef_name, const std::string &e_field_name)
    : hephaestus::ComplexMaxwellFormulation(
          magnetic_reluctivity_name, electric_conductivity_name,
          dielectric_permittivity_name, frequency_coef_name, e_field_name){};

} // namespace hephaestus

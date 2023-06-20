#include "complex_e_formulation.hpp"

namespace hephaestus {

ComplexEFormulation::ComplexEFormulation()
    : hephaestus::ComplexMaxwellFormulation() {
  h_curl_var_name = std::string("electric_field");

  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  zeta_coef_name = std::string("dielectric_permittivity");
};

void ComplexEFormulation::RegisterCoefficients() {
  hephaestus::DomainProperties &domain_properties =
      this->GetProblem()->domain_properties;

  freqCoef = dynamic_cast<mfem::ConstantCoefficient *>(
      domain_properties.scalar_property_map[frequency_coef_name]);
  if (freqCoef == NULL) {
    MFEM_ABORT("No frequency coefficient found. Frequency must be specified "
               "for frequency domain formulations.");
  }
  // define transformed
  domain_properties.scalar_property_map["_angular_frequency"] =
      new mfem::ConstantCoefficient(2.0 * M_PI * freqCoef->constant);
  domain_properties.scalar_property_map["_neg_angular_frequency"] =
      new mfem::ConstantCoefficient(-2.0 * M_PI * freqCoef->constant);
  domain_properties.scalar_property_map["_angular_frequency_sq"] =
      new mfem::ConstantCoefficient(pow(2.0 * M_PI * freqCoef->constant, 2));
  domain_properties.scalar_property_map["_neg_angular_frequency_sq"] =
      new mfem::ConstantCoefficient(-pow(2.0 * M_PI * freqCoef->constant, 2));

  if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
      0) {
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("magnetic_permeability")));
  }
  if (domain_properties.scalar_property_map.count(beta_coef_name) == 0) {
    domain_properties.scalar_property_map[beta_coef_name] =
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(beta_coef_name));
  }
  if (domain_properties.scalar_property_map.count(zeta_coef_name) == 0) {
    domain_properties.scalar_property_map[zeta_coef_name] =
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(zeta_coef_name));
  }

  domain_properties.scalar_property_map["maxwell_mass"] =
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map["_neg_angular_frequency_sq"],
          domain_properties.scalar_property_map[zeta_coef_name], prodFunc);

  domain_properties.scalar_property_map["maxwell_loss"] =
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map["_angular_frequency"],
          domain_properties.scalar_property_map[beta_coef_name], prodFunc);

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);
}

} // namespace hephaestus

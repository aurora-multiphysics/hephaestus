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
  hephaestus::Coefficients &domain_properties =
      this->GetProblem()->domain_properties;

  freqCoef = dynamic_cast<mfem::ConstantCoefficient *>(
      domain_properties.scalar_property_map.Get(frequency_coef_name));
  if (freqCoef == NULL) {
    MFEM_ABORT("No frequency coefficient found. Frequency must be specified "
               "for frequency domain formulations.");
  }
  // define transformed
  domain_properties.scalar_property_map.Register(
      "_angular_frequency",
      new mfem::ConstantCoefficient(2.0 * M_PI * freqCoef->constant), true);
  domain_properties.scalar_property_map.Register(
      "_neg_angular_frequency",
      new mfem::ConstantCoefficient(-2.0 * M_PI * freqCoef->constant), true);
  domain_properties.scalar_property_map.Register(
      "_angular_frequency_sq",
      new mfem::ConstantCoefficient(pow(2.0 * M_PI * freqCoef->constant, 2)),
      true);
  domain_properties.scalar_property_map.Register(
      "_neg_angular_frequency_sq",
      new mfem::ConstantCoefficient(-pow(2.0 * M_PI * freqCoef->constant, 2)),
      true);

  //   if (!domain_properties.scalar_property_map.Has("magnetic_permeability"))
  //   {
  //     domain_properties.scalar_property_map.Register(
  //         "magnetic_permeability",
  //         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //             std::string("magnetic_permeability"))),
  //         true);
  //   }
  //   if (!domain_properties.scalar_property_map.Has(beta_coef_name)) {
  //     domain_properties.scalar_property_map.Register(
  //         beta_coef_name,
  //         new mfem::PWCoefficient(
  //             domain_properties.getGlobalScalarProperty(beta_coef_name)),
  //         true);
  //   }
  //   if (!domain_properties.scalar_property_map.Has(zeta_coef_name)) {
  //     domain_properties.scalar_property_map.Register(
  //         zeta_coef_name,
  //         new mfem::PWCoefficient(
  //             domain_properties.getGlobalScalarProperty(zeta_coef_name)),
  //         true);
  //   }

  domain_properties.scalar_property_map.Register(
      "maxwell_mass",
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map.Get(
              "_neg_angular_frequency_sq"),
          domain_properties.scalar_property_map.Get(zeta_coef_name), prodFunc),
      true);

  domain_properties.scalar_property_map.Register(
      "maxwell_loss",
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map.Get("_angular_frequency"),
          domain_properties.scalar_property_map.Get(beta_coef_name), prodFunc),
      true);

  domain_properties.scalar_property_map.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map.Get("magnetic_permeability"),
          fracFunc),
      true);
}

} // namespace hephaestus

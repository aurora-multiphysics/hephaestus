#include "complex_a_form.hpp"

namespace hephaestus {

ComplexAFormulation::ComplexAFormulation() : ComplexMaxwellFormulation() {
  frequency_coef_name = std::string("frequency");
  h_curl_var_name = std::string("magnetic_vector_potential");

  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  zeta_coef_name = std::string("dielectric_permittivity");
};

void ComplexAFormulation::RegisterAuxSolvers(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::AuxSolvers &auxsolvers) {
  std::vector<std::string> aux_var_names;
  std::string b_field_name = "magnetic_flux_density";
  if (variables.Get(b_field_name + "_real") != NULL) {
    // if (myid_ == 0) {
    std::cout << b_field_name + "_real"
              << " found in variables: building auxvar " << std::endl;
    // }
    hephaestus::InputParameters b_field_aux_params;
    b_field_aux_params.SetParam("VariableName", h_curl_var_name + "_real");
    b_field_aux_params.SetParam("CurlVariableName", b_field_name + "_real");
    auxsolvers.Register("_magnetic_flux_density_re_aux",
                        new hephaestus::CurlAuxSolver(b_field_aux_params),
                        true);

    b_field_aux_params.SetParam("VariableName", h_curl_var_name + "_imag");
    b_field_aux_params.SetParam("CurlVariableName", b_field_name + "_imag");
    auxsolvers.Register("_magnetic_flux_density_im_aux",
                        new hephaestus::CurlAuxSolver(b_field_aux_params),
                        true);
  }
}

void ComplexAFormulation::RegisterCoefficients(
    hephaestus::DomainProperties &domain_properties) {

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

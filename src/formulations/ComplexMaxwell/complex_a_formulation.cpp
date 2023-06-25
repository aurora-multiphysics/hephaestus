#include "complex_a_formulation.hpp"

namespace hephaestus {

ComplexAFormulation::ComplexAFormulation() : ComplexMaxwellFormulation() {
  frequency_coef_name = std::string("frequency");
  h_curl_var_name = std::string("magnetic_vector_potential");

  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  zeta_coef_name = std::string("dielectric_permittivity");
};

void ComplexAFormulation::RegisterAuxSolvers() {
  hephaestus::DomainProperties &domain_properties =
      this->GetProblem()->domain_properties;
  hephaestus::GridFunctions &variables = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
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

  // E = -iωA
  std::string e_field_name = std::string("electric_field");
  if (variables.Get(e_field_name + "_real") != NULL) {
    // if (myid_ == 0) {
    std::cout << e_field_name + "_real"
              << " found in variables: building auxvar " << std::endl;
    // }
    hephaestus::InputParameters e_field_aux_params;
    e_field_aux_params.SetParam("CoefficientName",
                                std::string("_neg_angular_frequency"));
    e_field_aux_params.SetParam("InputVariableName", h_curl_var_name + "_real");
    e_field_aux_params.SetParam("ScaledVariableName", e_field_name + "_imag");
    auxsolvers.Register(
        "_electric_field_re_aux",
        new hephaestus::ScaledGridFunctionAuxSolver(e_field_aux_params), true);
    e_field_aux_params.SetParam("CoefficientName",
                                std::string("_angular_frequency"));
    e_field_aux_params.SetParam("InputVariableName", h_curl_var_name + "_imag");
    e_field_aux_params.SetParam("ScaledVariableName", e_field_name + "_real");
    auxsolvers.Register(
        "_electric_field_im_aux",
        new hephaestus::ScaledGridFunctionAuxSolver(e_field_aux_params), true);
  }
}

void ComplexAFormulation::RegisterCoefficients() {
  hephaestus::DomainProperties &domain_properties =
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

  if (!domain_properties.scalar_property_map.Has("magnetic_permeability")) {
    domain_properties.scalar_property_map.Register(
        "magnetic_permeability",
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("magnetic_permeability"))),
        true);
  }
  if (!domain_properties.scalar_property_map.Has(beta_coef_name)) {
    domain_properties.scalar_property_map.Register(
        beta_coef_name,
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(beta_coef_name)),
        true);
  }
  if (!domain_properties.scalar_property_map.Has(zeta_coef_name)) {
    domain_properties.scalar_property_map.Register(
        zeta_coef_name,
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(zeta_coef_name)),
        true);
  }

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

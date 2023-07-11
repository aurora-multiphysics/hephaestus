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
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  std::vector<std::string> aux_var_names;
  std::string b_field_name = "magnetic_flux_density";
  if (gridfunctions.Get(b_field_name + "_real") != NULL) {
    // if (myid_ == 0) {
    std::cout << b_field_name + "_real"
              << " found in gridfunctions: building auxvar " << std::endl;
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
  if (gridfunctions.Get(e_field_name + "_real") != NULL) {
    // if (myid_ == 0) {
    std::cout << e_field_name + "_real"
              << " found in gridfunctions: building auxvar " << std::endl;
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

} // namespace hephaestus

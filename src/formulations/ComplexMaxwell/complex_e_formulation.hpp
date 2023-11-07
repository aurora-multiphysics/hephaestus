#pragma once
#include "complex_maxwell_formulation.hpp"

namespace hephaestus {
/*
Formulation for solving:
∇×(ν∇×E꜀) + iωσE꜀ - ω²εE꜀ = -iωJ꜀ᵉ

via the weak form:
(ν∇×E꜀, ∇×E꜀') + (iωσE꜀, E꜀') - (ω²εE꜀, E꜀') - <(ν∇×E꜀)×n, E꜀'> = -(iωJ꜀ᵉ, E꜀')

where
E꜀ ∈ H(curl) is the complex electric field
J꜀ᵉ ∈ H(div) is the divergence-free source current density
ω is the angular frequency
ν is the reluctivity (reciprocal of magnetic permeability; 1/µ)
σ is the electric conductivity
ε is the electric permittivity

Dirichlet boundaries strongly constrain n×n×E꜀
Integrated boundaries weakly constrain (ν∇×E꜀)×n = n×dB꜀/dt
Robin boundaries weakly constrain (ν∇×E꜀)×n + γ(n×n×E꜀) = F

Divergence cleaning (such as via Helmholtz projection)
should be performed on J꜀ᵉ before use in this operator.
*/
class ComplexEFormulation : public hephaestus::ComplexMaxwellFormulation {
public:
  ComplexEFormulation(const std::string &magnetic_reluctivity_name,
                      const std::string &electric_conductivity_name,
                      const std::string &dielectric_permittivity_name,
                      const std::string &frequency_coef_name,
                      const std::string &e_field_complex_name,
                      const std::string &e_field_real_name,
                      const std::string &e_field_imag_name);

  // Enable auxiliary calculation of J ∈ H(div)
  virtual void
  registerCurrentDensityAux(const std::string &j_field_real_name,
                            const std::string &j_field_imag_name) override {
    //* Current density J = Jᵉ + σE
    //* Induced electric field, Jind = σE
    hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
    auxsolvers.Register(j_field_real_name,
                        new hephaestus::ScaledVectorGridFunctionAux(
                            _electric_field_real_name, j_field_real_name,
                            _electric_conductivity_name),
                        true);
  };

  // Enable auxiliary calculation of P ∈ L2
  virtual void registerJouleHeatingDensityAux(
      const std::string &p_field_name, const std::string &e_field_real_name,
      const std::string &e_field_imag_name,
      const std::string &conductivity_coef_name) override{

  };

protected:
  const std::string &_magnetic_reluctivity_name =
      hephaestus::ComplexMaxwellFormulation::_alpha_coef_name;
  const std::string &_electric_conductivity_name =
      hephaestus::ComplexMaxwellFormulation::_beta_coef_name;
  const std::string &_electric_permittivity_name =
      hephaestus::ComplexMaxwellFormulation::_zeta_coef_name;
  const std::string &_electric_field_complex_name =
      hephaestus::ComplexMaxwellFormulation::_h_curl_var_complex_name;
  const std::string &_electric_field_real_name =
      hephaestus::ComplexMaxwellFormulation::_h_curl_var_real_name;
  const std::string &_electric_field_imag_name =
      hephaestus::ComplexMaxwellFormulation::_h_curl_var_imag_name;
};
} // namespace hephaestus

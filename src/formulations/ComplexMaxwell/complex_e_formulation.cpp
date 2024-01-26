#include "complex_e_formulation.hpp"

namespace hephaestus
{

ComplexEFormulation::ComplexEFormulation(const std::string & magnetic_reluctivity_name,
                                         const std::string & electric_conductivity_name,
                                         const std::string & dielectric_permittivity_name,
                                         const std::string & frequency_coef_name,
                                         const std::string & e_field_complex_name,
                                         const std::string & e_field_real_name,
                                         const std::string & e_field_imag_name)
  : hephaestus::ComplexMaxwellFormulation(magnetic_reluctivity_name,
                                          electric_conductivity_name,
                                          dielectric_permittivity_name,
                                          frequency_coef_name,
                                          e_field_complex_name,
                                          e_field_real_name,
                                          e_field_imag_name)
{
}

// Enable auxiliary calculation of J ∈ H(div)
void
ComplexEFormulation::RegisterCurrentDensityAux(const std::string & j_field_real_name,
                                               const std::string & j_field_imag_name)
{
  //* Current density J = Jᵉ + σE
  //* Induced electric current, Jind = σE
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(
      j_field_real_name,
      std::make_shared<hephaestus::ScaledVectorGridFunctionAux>(
          _electric_field_real_name, j_field_real_name, _electric_conductivity_name));
  auxsolvers.Register(
      j_field_imag_name,
      std::make_shared<hephaestus::ScaledVectorGridFunctionAux>(
          _electric_field_imag_name, j_field_imag_name, _electric_conductivity_name));
};

void
ComplexEFormulation::RegisterMagneticFluxDensityAux(const std::string & b_field_real_name,
                                                    const std::string & b_field_imag_name)
{
  //* Magnetic flux density B = (i/ω) curl E
  //* (∇×E = -dB/dt = -iωB)
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(
      b_field_imag_name,
      std::make_shared<hephaestus::ScaledCurlVectorGridFunctionAux>(
          _electric_field_real_name, b_field_imag_name, "_inv_angular_frequency", 1.0));
  auxsolvers.Register(
      b_field_real_name,
      std::make_shared<hephaestus::ScaledCurlVectorGridFunctionAux>(
          _electric_field_imag_name, b_field_real_name, "_inv_angular_frequency", -1.0));
}

// Enable auxiliary calculation of P ∈ L2
void
ComplexEFormulation::RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                                    const std::string & e_field_real_name,
                                                    const std::string & e_field_imag_name,
                                                    const std::string & conductivity_coef_name)
{
  //* Time averaged Joule heating density = E.J
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(
      p_field_name,
      std::make_shared<hephaestus::VectorGridFunctionDotProductAux>(p_field_name,
                                                                    p_field_name,
                                                                    _electric_conductivity_name,
                                                                    _electric_field_real_name,
                                                                    _electric_field_real_name,
                                                                    _electric_field_imag_name,
                                                                    _electric_field_imag_name,
                                                                    true));
  auxsolvers.GetShared(p_field_name)->SetPriority(2);
}

} // namespace hephaestus

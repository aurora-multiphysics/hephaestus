#include "complex_a_formulation.hpp"

namespace hephaestus
{

ComplexAFormulation::ComplexAFormulation(const std::string & magnetic_reluctivity_name,
                                         const std::string & electric_conductivity_name,
                                         const std::string & dielectric_permittivity_name,
                                         const std::string & frequency_coef_name,
                                         const std::string & magnetic_vector_potential_complex_name,
                                         const std::string & magnetic_vector_potential_real_name,
                                         const std::string & magnetic_vector_potential_imag_name)
  : ComplexMaxwellFormulation(magnetic_reluctivity_name,
                              electric_conductivity_name,
                              dielectric_permittivity_name,
                              frequency_coef_name,
                              magnetic_vector_potential_complex_name,
                              magnetic_vector_potential_real_name,
                              magnetic_vector_potential_imag_name)
{
}

// Enable auxiliary calculation of J ∈ H(div)
void
ComplexAFormulation::RegisterCurrentDensityAux(const std::string & j_field_real_name,
                                               const std::string & j_field_imag_name)
{
  //* Current density J = Jᵉ + σE
  //* Induced current density Jind = σE = -iωσA
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(
      j_field_imag_name,
      std::make_shared<hephaestus::ScaledVectorGridFunctionAux>(
          _magnetic_vector_potential_real_name, j_field_imag_name, _loss_coef_name, -1.0));
  auxsolvers.Register(
      j_field_real_name,
      std::make_shared<hephaestus::ScaledVectorGridFunctionAux>(
          _magnetic_vector_potential_imag_name, j_field_real_name, _loss_coef_name, 1.0));
};

void
ComplexAFormulation::RegisterMagneticFluxDensityAux(const std::string & b_field_real_name,
                                                    const std::string & b_field_imag_name)
{
  //* Magnetic flux density B = curl A
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(b_field_real_name,
                      std::make_shared<hephaestus::CurlAuxSolver>(
                          _magnetic_vector_potential_real_name, b_field_real_name));
  auxsolvers.Register(b_field_imag_name,
                      std::make_shared<hephaestus::CurlAuxSolver>(
                          _magnetic_vector_potential_imag_name, b_field_imag_name));
}

void
ComplexAFormulation::RegisterElectricFieldAux(const std::string & e_field_real_name,
                                              const std::string & e_field_imag_name)
{
  //* Electric field E =-dA/dt=-iωA
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(
      e_field_imag_name,
      std::make_shared<hephaestus::ScaledVectorGridFunctionAux>(
          _magnetic_vector_potential_real_name, e_field_imag_name, "_angular_frequency", -1.0));
  auxsolvers.Register(
      e_field_real_name,
      std::make_shared<hephaestus::ScaledVectorGridFunctionAux>(
          _magnetic_vector_potential_imag_name, e_field_real_name, "_angular_frequency", 1.0));
}

// Enable auxiliary calculation of P ∈ L2
// Time averaged Joule heating density E.J
void
ComplexAFormulation::RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                                    const std::string & e_field_real_name,
                                                    const std::string & e_field_imag_name,
                                                    const std::string & conductivity_coef_name)
{
  //* Time averaged Joule heating density = E.J
  hephaestus::AuxSolvers & auxsolvers = GetProblem()->_postprocessors;
  auxsolvers.Register(p_field_name,
                      std::make_shared<hephaestus::VectorGridFunctionDotProductAux>(
                          p_field_name,
                          p_field_name,
                          _loss_coef_name,
                          _magnetic_vector_potential_real_name,
                          _magnetic_vector_potential_real_name,
                          _magnetic_vector_potential_imag_name,
                          _magnetic_vector_potential_imag_name,
                          true));

  auxsolvers.GetPtr(p_field_name)->SetPriority(2);
}

} // namespace hephaestus

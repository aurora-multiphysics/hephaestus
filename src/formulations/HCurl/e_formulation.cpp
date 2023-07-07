//* Solves:
//* ∇×(ν∇×E) + σdE/dt = -dJᵉ/dt
//*
//* in weak form
//* (ν∇×E, ∇×E') + (σdE/dt, E') + (dJᵉ/dt, E') - <(ν∇×E)×n, E'>  = 0

//* where:
//* reluctivity ν = 1/μ
//* electrical_conductivity σ=1/ρ
//* Electric Field E
//* Current density J = σE
//* Magnetic flux density, dB/dt = -∇×E
//* Magnetic field dH/dt = -ν∇×E

#include "e_formulation.hpp"

namespace hephaestus {

EFormulation::EFormulation() : HCurlFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  h_curl_var_name = std::string("electric_field");
}

void EFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &domain_properties =
      this->GetProblem()->domain_properties;
  //   if (!domain_properties.scalar_property_map.Has("magnetic_permeability"))
  //   {
  //     domain_properties.scalar_property_map.Register(
  //         "magnetic_permeability",
  //         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //             std::string("magnetic_permeability"))),
  //         true);
  //   }
  //   if
  //   (!domain_properties.scalar_property_map.Has("electrical_conductivity")) {
  //     domain_properties.scalar_property_map.Register(
  //         "electrical_conductivity",
  //         new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //             std::string("electrical_conductivity"))),
  //         true);
  //   }

  domain_properties.scalar_property_map.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map.Get("magnetic_permeability"),
          fracFunc),
      true);
}

} // namespace hephaestus

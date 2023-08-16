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

EFormulation::EFormulation(const std::string &magnetic_reluctivity_name,
                           const std::string &magnetic_permeability_name,
                           const std::string &electric_conductivity_name,
                           const std::string &e_field_name)
    : HCurlFormulation(magnetic_reluctivity_name, electric_conductivity_name,
                       e_field_name),
      _magnetic_permeability_name(magnetic_permeability_name) {}

void EFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has(_magnetic_permeability_name)) {
    MFEM_ABORT(_magnetic_permeability_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has(_electric_conductivity_name)) {
    MFEM_ABORT(_electric_conductivity_name + " coefficient not found.");
  }
  coefficients.scalars.Register(
      _magnetic_reluctivity_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get(_magnetic_permeability_name),
          fracFunc),
      true);
}
} // namespace hephaestus

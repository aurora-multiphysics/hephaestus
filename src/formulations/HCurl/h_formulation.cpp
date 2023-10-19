//* Solves:
//* ∇×(ρ∇×H) + μdH/dt = -dBᵉ/dt
//*
//* in weak form
//* (ρ∇×H, ∇×v) + (μdH/dt, v) + (dBᵉ/dt, v) - <(ρ∇×H)×n, v>  = 0

//* where:
//* magnetic permeability μ = 1/ν
//* electrical resistivity ρ=1/σ
//* Magnetic field, E = ρ∇×H
//* Magnetic flux density, B = Bᵉ + μH
//* Current density J = ∇×H

#include "h_formulation.hpp"

namespace hephaestus {

HFormulation::HFormulation(const std::string &electric_resistivity_name,
                           const std::string &electric_conductivity_name,
                           const std::string &magnetic_permeability_name,
                           const std::string &h_field_name)
    : HCurlFormulation(electric_resistivity_name, magnetic_permeability_name,
                       h_field_name),
      _electric_conductivity_name(electric_conductivity_name) {}

void HFormulation::RegisterAuxSolvers() {
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  std::vector<std::string> aux_var_names;
  std::string j_field_name = "current_density";
  if (gridfunctions.Get(j_field_name) != NULL) {
    auxsolvers.Register(
        "_current_density_aux",
        new hephaestus::CurlAuxSolver(_h_curl_var_name, j_field_name), true);
  }
}

void HFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has(_electric_conductivity_name)) {
    MFEM_ABORT(_electric_conductivity_name + " coefficient not found.");
  }
  coefficients.scalars.Register(
      _electric_resistivity_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get(_electric_conductivity_name),
          fracFunc),
      true);
  HCurlFormulation::RegisterCoefficients();
}

} // namespace hephaestus

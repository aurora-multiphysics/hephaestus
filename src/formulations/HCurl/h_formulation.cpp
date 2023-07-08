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

HFormulation::HFormulation() : HCurlFormulation() {
  alpha_coef_name = std::string("electrical_resistivity");
  beta_coef_name = std::string("magnetic_permeability");
  h_curl_var_name = std::string("magnetic_field");
}

void HFormulation::RegisterAuxSolvers() {
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  std::vector<std::string> aux_var_names;
  std::string j_field_name = "current_density";
  if (gridfunctions.Get(j_field_name) != NULL) {
    hephaestus::InputParameters j_field_aux_params;
    j_field_aux_params.SetParam("VariableName", std::string("magnetic_field"));
    j_field_aux_params.SetParam("CurlVariableName", j_field_name);
    auxsolvers.Register("_magnetic_flux_density_aux",
                        new hephaestus::CurlAuxSolver(j_field_aux_params),
                        true);
  }
}

void HFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has("magnetic_permeability")) {
    MFEM_ABORT("magnetic_permeability coefficient not found.");
  }
  if (!coefficients.scalars.Has("electrical_conductivity")) {
    MFEM_ABORT("electrical_conductivity coefficient not found.");
  }

  coefficients.scalars.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get("electrical_conductivity"),
          fracFunc),
      true);
}

} // namespace hephaestus

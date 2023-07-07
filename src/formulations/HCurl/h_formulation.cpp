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
  hephaestus::GridFunctions &variables = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  std::vector<std::string> aux_var_names;
  std::string j_field_name = "current_density";
  if (variables.Get(j_field_name) != NULL) {
    hephaestus::InputParameters j_field_aux_params;
    j_field_aux_params.SetParam("VariableName", std::string("magnetic_field"));
    j_field_aux_params.SetParam("CurlVariableName", j_field_name);
    auxsolvers.Register("_magnetic_flux_density_aux",
                        new hephaestus::CurlAuxSolver(j_field_aux_params),
                        true);
  }
}

void HFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &domain_properties =
      this->GetProblem()->domain_properties;
  // if (!domain_properties.scalar_property_map.Has("magnetic_permeability")) {
  //   domain_properties.scalar_property_map.Register(
  //       "magnetic_permeability",
  //       new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //           std::string("magnetic_permeability"))),
  //       true);
  // }
  // if (!domain_properties.scalar_property_map.Has("electrical_conductivity"))
  // {
  //   domain_properties.scalar_property_map.Register(
  //       "electrical_conductivity",
  //       new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //           std::string("electrical_conductivity"))),
  //       true);
  // }

  domain_properties.scalar_property_map.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map.Get("electrical_conductivity"),
          fracFunc),
      true);
}

} // namespace hephaestus

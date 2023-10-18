//* Solves:
//* ∇×(ν∇×A) = Jᵉ
//*
//* in weak form
//* (ν∇×A, ∇×A') - (Jᵉ, A') - <(ν∇×A)×n, A'>  = 0

//* where:
//* reluctivity ν = 1/μ
//* Magnetic vector potential A
//* Electric field, E = ρJᵉ
//* Magnetic flux density, B = ∇×A
//* Magnetic field H = ν∇×A
//* Current density J = Jᵉ

#include "magnetostatic_formulation.hpp"

namespace hephaestus {

MagnetostaticFormulation::MagnetostaticFormulation(
    const std::string &magnetic_reluctivity_name,
    const std::string &magnetic_permeability_name,
    const std::string &magnetic_vector_potential_name)
    : StaticsFormulation(magnetic_reluctivity_name,
                         magnetic_vector_potential_name),
      _magnetic_permeability_name(magnetic_permeability_name) {}

void MagnetostaticFormulation::RegisterAuxSolvers() {
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  std::vector<std::string> aux_var_names;
  std::string b_field_name = "magnetic_flux_density";
  if (gridfunctions.Get(b_field_name) != NULL) {
    hephaestus::InputParameters b_field_aux_params;
    b_field_aux_params.SetParam("VariableName", _h_curl_var_name);
    b_field_aux_params.SetParam("CurlVariableName", b_field_name);
    auxsolvers.Register("_magnetic_flux_density_aux",
                        new hephaestus::CurlAuxSolver(b_field_aux_params),
                        true);
  }
}

void MagnetostaticFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has(_magnetic_permeability_name)) {
    MFEM_ABORT(_magnetic_permeability_name + " coefficient not found.");
  }
  coefficients.scalars.Register(
      _magnetic_reluctivity_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get(_magnetic_permeability_name),
          fracFunc),
      true);
}
} // namespace hephaestus

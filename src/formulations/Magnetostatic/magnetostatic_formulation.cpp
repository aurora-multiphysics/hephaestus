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

void MagnetostaticFormulation::registerMagneticFluxDensityAux(
    const std::string &b_field_name) {
  //* Magnetic flux density, B = ∇×A
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      b_field_name,
      new hephaestus::CurlAuxSolver(_h_curl_var_name, b_field_name), true);
}

void MagnetostaticFormulation::registerMagneticFieldAux(
    const std::string &h_field_name) {
  //* Magnetic field H = ν∇×A
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      h_field_name,
      new hephaestus::ScaledCurlVectorGridFunctionAux(
          _h_curl_var_name, h_field_name, _magnetic_reluctivity_name),
      true);
}

void MagnetostaticFormulation::registerLorentzForceDensityAux(
    const std::string &f_field_name, const std::string &b_field_name,
    const std::string &j_field_name) {
  //* Lorentz force density = J x B
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      f_field_name,
      new hephaestus::VectorGridFunctionCrossProductAux(
          f_field_name, f_field_name, j_field_name, b_field_name),
      true);
  auxsolvers.Get(f_field_name)->SetPriority(2);
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

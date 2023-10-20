//* Solves:
//* ∇×(ρ∇×H) + μdH/dt = -dBᵉ/dt
//*
//* in weak form
//* (ρ∇×H, ∇×v) + (μdH/dt, v) + (dBᵉ/dt, v) - <(ρ∇×H)×n, v>  = 0

//* where:
//* magnetic permeability μ = 1/ν
//* electrical resistivity ρ=1/σ
//* Electric field, E = ρ∇×H
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

void HFormulation::registerCurrentDensityAux(const std::string &j_field_name) {
  //* Current density J = ∇×H
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      j_field_name,
      new hephaestus::CurlAuxSolver(_h_curl_var_name, j_field_name), true);
}

void HFormulation::registerMagneticFluxDensityAux(
    const std::string &b_field_name) {
  //* Magnetic flux density, B = Bᵉ + μH
  //* Induced flux density, B = μH
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      b_field_name,
      new hephaestus::ScaledVectorGridFunctionAux(
          _h_curl_var_name, b_field_name, _magnetic_permeability_name),
      true);
}

void HFormulation::registerElectricFieldAux(const std::string &e_field_name) {
  //* Electric field, E = ρ∇×H
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      e_field_name,
      new hephaestus::ScaledCurlVectorGridFunctionAux(
          _h_curl_var_name, e_field_name, _electric_resistivity_name),
      true);
}

void HFormulation::registerMagneticFieldAux(const std::string &h_field_name) {
  //* Magnetic field H is a state variable; no additional calculation needed
}

void HFormulation::registerLorentzForceDensityAux(
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

void HFormulation::registerJouleHeatingDensityAux(
    const std::string &p_field_name, const std::string &e_field_name,
    const std::string &conductivity_coef_name) {
  //* Joule heating density = E.J
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(p_field_name,
                      new hephaestus::VectorGridFunctionDotProductAux(
                          p_field_name, p_field_name, e_field_name,
                          e_field_name, _electric_conductivity_name),
                      true);
  auxsolvers.Get(p_field_name)->SetPriority(2);
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

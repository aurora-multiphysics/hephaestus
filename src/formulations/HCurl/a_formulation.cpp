//* Solves:
//* ∇×(ν∇×A) + σdA/dt = Jᵉ
//*
//* in weak form
//* (ν∇×A, ∇×A') + (σdA/dt, A') - (Jᵉ, A') - <(ν∇×A)×n, A'>  = 0

//* where:
//* reluctivity ν = 1/μ
//* electrical_conductivity σ=1/ρ
//* Magnetic vector potential A
//* Electric field, E = ρJᵉ -dA/dt
//* Magnetic flux density, B = ∇×A
//* Magnetic field H = ν∇×A
//* Current density J = Jᵉ -σdA/dt

#include "a_formulation.hpp"

namespace hephaestus {

AFormulation::AFormulation(const std::string &magnetic_reluctivity_name,
                           const std::string &magnetic_permeability_name,
                           const std::string &electric_conductivity_name,
                           const std::string &magnetic_vector_potential_name)
    : HCurlFormulation(magnetic_reluctivity_name, electric_conductivity_name,
                       magnetic_vector_potential_name),
      _magnetic_permeability_name(magnetic_permeability_name) {}

void AFormulation::RegisterAuxSolvers() {
  // registerCurrentDensityAux("current_density");
  // registerMagneticFluxDensityAux("magnetic_flux_density");
  // registerElectricFieldAux("electric_field");
  // registerMagneticFieldAux("magnetic_field");
}

void AFormulation::registerCurrentDensityAux(const std::string &j_field_name) {
  //* Current density J = Jᵉ -σdA/dt
  //* Induced electric field, Jind = -σdA/dt
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(j_field_name,
                      new hephaestus::ScaledVectorGridFunctionAux(
                          GetTimeDerivativeName(_h_curl_var_name), j_field_name,
                          _electric_conductivity_name, -1.0),
                      true);
}

void AFormulation::registerMagneticFluxDensityAux(
    const std::string &b_field_name) {
  //* Magnetic flux density, B = ∇×A
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      b_field_name,
      new hephaestus::CurlAuxSolver(_h_curl_var_name, b_field_name), true);
}

void AFormulation::registerElectricFieldAux(const std::string &e_field_name) {
  //* Total electric field, E = ρJᵉ -dA/dt
  //* Induced electric field, Eind = -dA/dt
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      e_field_name,
      new hephaestus::ScaledVectorGridFunctionAux(
          GetTimeDerivativeName(_h_curl_var_name), e_field_name, "_one", -1.0),
      true);
}

void AFormulation::registerMagneticFieldAux(const std::string &h_field_name) {
  //* Magnetic field H = ν∇×A
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      h_field_name,
      new hephaestus::ScaledCurlVectorGridFunctionAux(
          _h_curl_var_name, h_field_name, _magnetic_reluctivity_name),
      true);
}

void AFormulation::registerLorentzForceDensityAux(
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

void AFormulation::registerJouleHeatingDensityAux(
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

void AFormulation::RegisterCoefficients() {
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

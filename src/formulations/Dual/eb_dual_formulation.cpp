#include "eb_dual_formulation.hpp"

namespace hephaestus {

EBDualFormulation::EBDualFormulation(
    const std::string &magnetic_reluctivity_name,
    const std::string &magnetic_permeability_name,
    const std::string &electric_conductivity_name,
    const std::string &e_field_name, const std::string &b_field_name)
    : DualFormulation(magnetic_reluctivity_name, electric_conductivity_name,
                      e_field_name, b_field_name),
      _magnetic_permeability_name(magnetic_permeability_name) {}

void EBDualFormulation::registerCurrentDensityAux(
    const std::string &j_field_name) {
  //* Current density J = Jᵉ + σE
  //* Induced electric field, Jind = σE
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(
      j_field_name,
      new hephaestus::ScaledVectorGridFunctionAux(
          _h_curl_var_name, j_field_name, _electric_conductivity_name),
      true);
}

void EBDualFormulation::registerLorentzForceDensityAux(
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

void EBDualFormulation::registerJouleHeatingDensityAux(
    const std::string &p_field_name, const std::string &e_field_name,
    const std::string &conductivity_coef_name) {
  //* Joule heating density = E.J
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  auxsolvers.Register(p_field_name,
                      new hephaestus::VectorGridFunctionDotProductAux(
                          p_field_name, p_field_name,
                          _electric_conductivity_name, e_field_name,
                          e_field_name),
                      true);
  auxsolvers.Get(p_field_name)->SetPriority(2);
}

void EBDualFormulation::RegisterCoefficients() {
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
  DualFormulation::RegisterCoefficients();
}

} // namespace hephaestus

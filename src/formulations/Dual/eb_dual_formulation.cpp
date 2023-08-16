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

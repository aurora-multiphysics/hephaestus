#include "hj_dual_formulation.hpp"

namespace hephaestus {

HJDualFormulation::HJDualFormulation(
    const std::string &electric_resistivity_name,
    const std::string &electric_conductivity_name,
    const std::string &magnetic_permeability_name,
    const std::string &h_field_name, const std::string &j_field_name)
    : DualFormulation(electric_resistivity_name, magnetic_permeability_name,
                      h_field_name, j_field_name),
      _electric_conductivity_name(electric_conductivity_name) {}

void HJDualFormulation::RegisterCoefficients() {
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
  DualFormulation::RegisterCoefficients();
}

} // namespace hephaestus

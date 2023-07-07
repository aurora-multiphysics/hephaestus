#include "coefficient_aux.hpp"

namespace hephaestus {

CoefficientAuxSolver::CoefficientAuxSolver(
    const hephaestus::InputParameters &params)
    : AuxSolver(), var_name(params.GetParam<std::string>("VariableName")),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void CoefficientAuxSolver::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::Coefficients &domain_properties) {
  gf = variables.Get(var_name);
  if (gf == NULL) {
    MFEM_ABORT("GridFunction "
               << var_name
               << " not found when initializing CoefficientAuxSolver");
  }
  coeff = domain_properties.scalar_property_map.Get(coef_name);
  if (gf == NULL) {
    MFEM_ABORT("Coefficient "
               << coef_name
               << " not found when initializing CoefficientAuxSolver");
  }
}

void CoefficientAuxSolver::Solve(double t) {
  coeff->SetTime(t);
  gf->ProjectCoefficient(*coeff);
}

} // namespace hephaestus

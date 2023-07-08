#include "vector_coefficient_aux.hpp"

namespace hephaestus {

VectorCoefficientAuxSolver::VectorCoefficientAuxSolver(
    const hephaestus::InputParameters &params)
    : AuxSolver(), var_name(params.GetParam<std::string>("VariableName")),
      vec_coef_name(params.GetParam<std::string>("VectorCoefficientName")) {}

void VectorCoefficientAuxSolver::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::Coefficients &coefficients) {
  gf = variables.Get(var_name);
  if (gf == NULL) {
    MFEM_ABORT("GridFunction "
               << var_name
               << " not found when initializing CoefficientAuxSolver");
  }
  vec_coeff = coefficients.vectors.Get(vec_coef_name);
  if (vec_coeff == NULL) {
    MFEM_ABORT("VectorCoefficient "
               << vec_coef_name
               << " not found when initializing VectorCoefficientAuxSolver");
  }
}

void VectorCoefficientAuxSolver::Solve(double t) {
  vec_coeff->SetTime(t);
  gf->ProjectCoefficient(*vec_coeff);
}

} // namespace hephaestus

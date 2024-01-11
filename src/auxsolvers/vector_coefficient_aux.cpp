#include "vector_coefficient_aux.hpp"

namespace hephaestus {

VectorCoefficientAux::VectorCoefficientAux(const std::string &gf_name,
                                           const std::string &vec_coef_name)
    : AuxSolver(), _gf_name(gf_name), _vec_coef_name(vec_coef_name) {}

void VectorCoefficientAux::Init(const hephaestus::GridFunctions &gridfunctions,
                                hephaestus::Coefficients &coefficients) {
  gf = gridfunctions.Get(_gf_name);
  if (gf == nullptr) {
    MFEM_ABORT("GridFunction "
               << _gf_name
               << " not found when initializing VectorCoefficientAux");
  }
  vec_coef = coefficients.vectors.Get(_vec_coef_name);
  if (vec_coef == nullptr) {
    MFEM_ABORT("VectorCoefficient "
               << _vec_coef_name
               << " not found when initializing VectorCoefficientAux");
  }
}

void VectorCoefficientAux::Solve(double t) {
  vec_coef->SetTime(t);
  gf->ProjectCoefficient(*vec_coef);
}

} // namespace hephaestus

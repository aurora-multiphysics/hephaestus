#include "coefficient_aux.hpp"

namespace hephaestus {

CoefficientAux::CoefficientAux(const std::string &gf_name,
                               const std::string &coef_name)
    : AuxSolver(), _gf_name(gf_name), _coef_name(coef_name) {}

void CoefficientAux::Init(const hephaestus::GridFunctions &gridfunctions,
                          hephaestus::Coefficients &coefficients) {
  gf = gridfunctions.Get(_gf_name);
  if (gf == NULL) {
    MFEM_ABORT("GridFunction "
               << _gf_name << " not found when initializing CoefficientAux");
  }
  coef = coefficients.scalars.Get(_coef_name);
  if (coef == NULL) {
    MFEM_ABORT("Coefficient " << _coef_name
                              << " not found when initializing CoefficientAux");
  }
}

void CoefficientAux::Solve(double t) { gf->ProjectCoefficient(*coef); }

} // namespace hephaestus

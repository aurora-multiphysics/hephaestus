#include "coefficient_aux.hpp"

namespace hephaestus
{

CoefficientAux::CoefficientAux(const std::string & gf_name, const std::string & coef_name)
  : AuxSolver(), _gf_name(gf_name), _coef_name(coef_name)
{
}

void
CoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  if (_gf == NULL)
  {
    MFEM_ABORT("GridFunction " << _gf_name << " not found when initializing CoefficientAux");
  }
  _coef = coefficients.scalars.Get(_coef_name);
  if (_coef == NULL)
  {
    MFEM_ABORT("Coefficient " << _coef_name << " not found when initializing CoefficientAux");
  }
}

void
CoefficientAux::Solve(double t)
{
  _gf->ProjectCoefficient(*_coef);
}

} // namespace hephaestus

#include "coefficient_aux.hpp"

#include <utility>

namespace hephaestus
{

CoefficientAux::CoefficientAux(std::string  gf_name, std::string  coef_name)
  :  _gf_name(std::move(gf_name)), _coef_name(std::move(coef_name))
{
}

void
CoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  if (_gf == nullptr)
  {
    MFEM_ABORT("GridFunction " << _gf_name << " not found when initializing CoefficientAux");
  }
  _coef = coefficients.scalars.Get(_coef_name);
  if (_coef == nullptr)
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

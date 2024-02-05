#include "coefficient_aux.hpp"

#include <utility>

namespace hephaestus
{

CoefficientAux::CoefficientAux(std::string gf_name, std::string coef_name)
  : _gf_name(std::move(gf_name)), _coef_name(std::move(coef_name))
{
}

void
CoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  // NB: ensure pointers are not NULL with "nullable = false".
  _gf = gridfunctions.GetPtr(_gf_name, false);
  _coef = coefficients._scalars.GetPtr(_coef_name, false);
}

void
CoefficientAux::Solve(double t)
{
  _gf->ProjectCoefficient(*_coef);
}

} // namespace hephaestus

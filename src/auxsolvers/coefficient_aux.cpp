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
  _gf = gridfunctions.Get(_gf_name);
  _coef = coefficients._scalars.Get(_coef_name);
}

void
CoefficientAux::Solve(double t)
{
  _gf->ProjectCoefficient(*_coef);
}

} // namespace hephaestus

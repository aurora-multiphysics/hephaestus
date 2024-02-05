#include "vector_coefficient_aux.hpp"

#include <utility>

namespace hephaestus
{

VectorCoefficientAux::VectorCoefficientAux(std::string gf_name, std::string vec_coef_name)
  : _gf_name(std::move(gf_name)), _vec_coef_name(std::move(vec_coef_name))
{
}

void
VectorCoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                           hephaestus::Coefficients & coefficients)
{
  // NB: set "nullable = false" to ensure pointers valid.
  _gf = gridfunctions.GetPtr(_gf_name, false);
  _vec_coef = coefficients._vectors.GetPtr(_vec_coef_name, false);
}

void
VectorCoefficientAux::Solve(double t)
{
  _vec_coef->SetTime(t);
  _gf->ProjectCoefficient(*_vec_coef);
}

} // namespace hephaestus

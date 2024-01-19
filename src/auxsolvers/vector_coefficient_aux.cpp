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
  gf = gridfunctions.Get(_gf_name);
  if (gf == nullptr)
  {
    MFEM_ABORT("GridFunction " << _gf_name << " not found when initializing VectorCoefficientAux");
  }
  vec_coef = coefficients.vectors.Get(_vec_coef_name);
  if (vec_coef == nullptr)
  {
    MFEM_ABORT("VectorCoefficient " << _vec_coef_name
                                    << " not found when initializing VectorCoefficientAux");
  }
}

void
VectorCoefficientAux::Solve(double t)
{
  vec_coef->SetTime(t);
  gf->ProjectCoefficient(*vec_coef);
}

} // namespace hephaestus

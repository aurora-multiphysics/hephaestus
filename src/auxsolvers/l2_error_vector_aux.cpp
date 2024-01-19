#include "l2_error_vector_aux.hpp"

namespace hephaestus
{

L2ErrorVectorPostprocessor::L2ErrorVectorPostprocessor(const hephaestus::InputParameters & params)
  : 
    var_name(params.GetParam<std::string>("VariableName")),
    vec_coef_name(params.GetParam<std::string>("VectorCoefficientName"))
{
}

void
L2ErrorVectorPostprocessor::Init(const hephaestus::GridFunctions & gridfunctions,
                                 hephaestus::Coefficients & coefficients)
{
  gf = gridfunctions.Get(var_name);
  vec_coeff = coefficients.vectors.Get(vec_coef_name);
}

void
L2ErrorVectorPostprocessor::Solve(double t)
{
  double l2_err = gf->ComputeL2Error(*vec_coeff);
  HYPRE_BigInt ndof = gf->ParFESpace()->GlobalTrueVSize();

  times.Append(t);
  l2_errs.Append(l2_err);
  ndofs.Append(ndof);
}

} // namespace hephaestus

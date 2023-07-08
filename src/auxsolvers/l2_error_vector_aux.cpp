#include "postprocessors.hpp"

namespace hephaestus {

L2ErrorVectorPostprocessor::L2ErrorVectorPostprocessor(
    const hephaestus::InputParameters &params)
    : AuxSolver(), var_name(params.GetParam<std::string>("VariableName")),
      vec_coef_name(params.GetParam<std::string>("VectorCoefficientName")) {}

void L2ErrorVectorPostprocessor::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::Coefficients &coefficients) {
  gf = variables.Get(var_name);
  vec_coeff = coefficients.vectors.Get(vec_coef_name);
}

void L2ErrorVectorPostprocessor::Solve(double t) {
  double l2_err = gf->ComputeL2Error(*vec_coeff);
  HYPRE_BigInt ndof = gf->ParFESpace()->GlobalTrueVSize();

  times.Append(t);
  l2_errs.Append(l2_err);
  ndofs.Append(ndof);
}

} // namespace hephaestus

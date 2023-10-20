#include "vector_gridfunction_cross_product_aux.hpp"

namespace hephaestus {

void cross_product(mfem::Vector &va, mfem::Vector &vb, mfem::Vector &V) {
  V.SetSize(3);
  V[0] = va[1] * vb[2] - va[2] * vb[1];
  V[1] = va[2] * vb[0] - va[0] * vb[2];
  V[2] = va[0] * vb[1] - va[1] * vb[0];
}

void VectorGridFunctionCrossProductCoefficient::Eval(
    mfem::Vector &uxv, mfem::ElementTransformation &T,
    const mfem::IntegrationPoint &ip) {
  mfem::Vector u_vec, v_vec;
  _u_gf.GetVectorValue(T, ip, u_vec);
  _v_gf.GetVectorValue(T, ip, v_vec);
  hephaestus::cross_product(u_vec, v_vec, uxv);
}

VectorGridFunctionCrossProductAux::VectorGridFunctionCrossProductAux(
    const std::string &cross_product_gf_name,
    const std::string &cross_product_coef_name, const std::string &u_gf_name,
    const std::string &v_gf_name)
    : VectorCoefficientAux(cross_product_gf_name, cross_product_coef_name),
      _u_gf_name(u_gf_name), _v_gf_name(v_gf_name), _u_gf(nullptr),
      _v_gf(nullptr) {}

void VectorGridFunctionCrossProductAux::Init(
    const hephaestus::GridFunctions &gridfunctions,
    hephaestus::Coefficients &coefficients) {
  _u_gf = gridfunctions.Get(_u_gf_name);
  if (_u_gf == NULL) {
    MFEM_ABORT(
        "GridFunction "
        << _u_gf_name
        << " not found when initializing VectorGridFunctionCrossProductAux");
  }
  _v_gf = gridfunctions.Get(_v_gf_name);
  if (_v_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << _v_gf_name
               << " not found when initializing ScaledVectorGridFunctionAux");
  }

  coefficients.vectors.Register(
      _vec_coef_name,
      new hephaestus::VectorGridFunctionCrossProductCoefficient(*_u_gf, *_v_gf),
      true);

  VectorCoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus

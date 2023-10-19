#include "lorentz_force_aux.hpp"

namespace hephaestus {

void cross_product(mfem::Vector &va, mfem::Vector &vb, mfem::Vector &V) {
  V.SetSize(3);
  V[0] = va[1] * vb[2] - va[2] * vb[1];
  V[1] = va[2] * vb[0] - va[0] * vb[2];
  V[2] = va[0] * vb[1] - va[1] * vb[0];
}

void LorentzForceDensityAuxCoefficient::Eval(mfem::Vector &V,
                                             mfem::ElementTransformation &T,
                                             const mfem::IntegrationPoint &ip) {
  mfem::Vector B, J;
  B_gf.GetVectorValue(T, ip, B);
  J_gf.GetVectorValue(T, ip, J);
  hephaestus::cross_product(J, B, V);
}

LorentzForceDensityAux::LorentzForceDensityAux(
    const std::string &lorentz_force_density_gf_name,
    const std::string &lorentz_force_density_coef_name,
    const std::string &magnetic_flux_density_name,
    const std::string &current_density_name)
    : VectorCoefficientAux(lorentz_force_density_gf_name,
                           lorentz_force_density_coef_name),
      _magnetic_flux_density_name(magnetic_flux_density_name),
      _current_density_name(current_density_name), j_gf(nullptr),
      b_gf(nullptr) {}

void LorentzForceDensityAux::Init(
    const hephaestus::GridFunctions &gridfunctions,
    hephaestus::Coefficients &coefficients) {
  j_gf = gridfunctions.Get(_current_density_name);
  if (j_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << _current_density_name
               << " not found when initializing LorentzForceDensityAux");
  }
  b_gf = gridfunctions.Get(_magnetic_flux_density_name);
  if (b_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << _magnetic_flux_density_name
               << " not found when initializing ScaledVectorGridFunctionAux");
  }

  coefficients.vectors.Register(
      _vec_coef_name,
      new hephaestus::LorentzForceDensityAuxCoefficient(*j_gf, *b_gf), true);

  VectorCoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus

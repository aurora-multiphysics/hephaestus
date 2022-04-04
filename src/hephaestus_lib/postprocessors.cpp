#include "postprocessors.hpp"

namespace hephaestus {

void cross_product(mfem::Vector &va, mfem::Vector &vb, mfem::Vector &V) {
  V.SetSize(3);
  V[0] = va[1] * vb[2] - va[2] * vb[1];
  V[1] = va[2] * vb[0] - va[0] * vb[2];
  V[2] = va[0] * vb[1] - va[1] * vb[0];
}

double JouleHeatingCoefficient::Eval(mfem::ElementTransformation &T,
                                     const mfem::IntegrationPoint &ip) {
  mfem::Vector E;
  double thisSigma;
  E_gf.GetVectorValue(T, ip, E);
  thisSigma = sigma.Eval(T, ip);
  return thisSigma * (E * E);
}

void LorentzForceVectorCoefficient::Eval(mfem::Vector &V,
                                         mfem::ElementTransformation &T,
                                         const mfem::IntegrationPoint &ip) {
  mfem::Vector B, J;

  B_gf.GetVectorValue(T, ip, B);
  J_gf.GetVectorValue(T, ip, J);

  hephaestus::cross_product(J, B, V);
}

} // namespace hephaestus
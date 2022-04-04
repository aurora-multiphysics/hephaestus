#include "postprocessors.hpp"

namespace hephaestus {

double JouleHeatingCoefficient::Eval(mfem::ElementTransformation &T,
                                     const mfem::IntegrationPoint &ip) {
  mfem::Vector E;
  double thisSigma;
  E_gf.GetVectorValue(T, ip, E);
  thisSigma = sigma.Eval(T, ip);
  return thisSigma * (E * E);
}

} // namespace hephaestus
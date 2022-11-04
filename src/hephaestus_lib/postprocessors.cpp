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

void Postprocessors::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  for (const auto &[name, postprocessor] : GetMap()) {
    postprocessor->Init(variables, domain_properties);
  }
}
void Postprocessors::Update(double t) {
  for (const auto &[name, postprocessor] : GetMap()) {
    postprocessor->Update(t);
  }
}

L2ErrorVectorPostprocessor::L2ErrorVectorPostprocessor(
    const hephaestus::InputParameters &params)
    : Postprocessor(params),
      var_name(params.GetParam<std::string>("VariableName")),
      vec_coef_name(params.GetParam<std::string>("VectorCoefficientName")) {}

void L2ErrorVectorPostprocessor::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  gf = variables.Get(var_name);
  vec_coeff = domain_properties.vector_property_map[vec_coef_name];
}

void L2ErrorVectorPostprocessor::Update(double t) {
  double l2_err = gf->ComputeL2Error(*vec_coeff);
  HYPRE_BigInt ndof = gf->ParFESpace()->GlobalTrueVSize();

  times.Append(t);
  l2_errs.Append(l2_err);
  ndofs.Append(ndof);
}

} // namespace hephaestus

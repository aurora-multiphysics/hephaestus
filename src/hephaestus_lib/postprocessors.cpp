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
// Constructor:
//  - name one or more variable names used in coef eval.
//  - provide name of coef
// Init: pass variables and domain properties, and assign ptrs.
//   - Create coef if not present.
// Eval: return coeff value at point

// Inherit from mfem::Coefficient or PostProcessor?
// Logic:
// Just need a PostProcessor that names coupled variable in constructor,
// adds its coef to domainproperties at init,
// and updates at update
//
// input: coupled var name, coef name, function

CoupledCoefficient::CoupledCoefficient(
    const hephaestus::InputParameters &params)
    : Postprocessor(params),
      coupled_var_name(params.GetParam<std::string>("CoupledVariableName")) {}

void CoupledCoefficient::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  if (variables.Has(coupled_var_name)) {
    gf = variables.Get(coupled_var_name);
  } else {
    const std::string error_message = coupled_var_name +
                                      " not found in variables when "
                                      "creating CoupledCoefficient\n";
    mfem::mfem_error(error_message.c_str());
  }
}

double CoupledCoefficient::Eval(mfem::ElementTransformation &T,
                                const mfem::IntegrationPoint &ip) {
  return gf->GetValue(T, ip);
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

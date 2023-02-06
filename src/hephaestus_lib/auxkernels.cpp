#include "auxkernels.hpp"

namespace hephaestus {

void AuxKernels::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  for (const auto &[name, auxkernel] : GetMap()) {
    auxkernel->Init(variables, domain_properties);
  }
}
void AuxKernels::Solve(double t) {
  for (const auto &[name, auxkernel] : GetMap()) {
    auxkernel->Solve(t);
  }
}

CurlAuxKernel::CurlAuxKernel(const hephaestus::InputParameters &params)
    : AuxKernel(params), var_name(params.GetParam<std::string>("VariableName")),
      curl_var_name(params.GetParam<std::string>("CurlVariableName")) {}

void CurlAuxKernel::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  u_ = variables.Get(var_name);
  curl_u_ = variables.Get(curl_var_name);

  curl = new mfem::ParDiscreteLinearOperator(u_->ParFESpace(),
                                             curl_u_->ParFESpace());
  curl->AddDomainInterpolator(new mfem::CurlInterpolator());
  curl->Assemble();
}

void CurlAuxKernel::Solve(double t) { curl->Mult(*u_, *curl_u_); }

CoefficientAuxKernel::CoefficientAuxKernel(
    const hephaestus::InputParameters &params)
    : AuxKernel(params), var_name(params.GetParam<std::string>("VariableName")),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void CoefficientAuxKernel::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  gf = variables.Get(var_name);
  coeff = domain_properties.scalar_property_map[coef_name];
}

void CoefficientAuxKernel::Solve(double t) {
  coeff->SetTime(t);
  gf->ProjectCoefficient(*coeff);
}

VectorCoefficientAuxKernel::VectorCoefficientAuxKernel(
    const hephaestus::InputParameters &params)
    : AuxKernel(params), var_name(params.GetParam<std::string>("VariableName")),
      vec_coef_name(params.GetParam<std::string>("VectorCoefficientName")) {}

void VectorCoefficientAuxKernel::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  gf = variables.Get(var_name);
  vec_coeff = domain_properties.vector_property_map[vec_coef_name];
}

void VectorCoefficientAuxKernel::Solve(double t) {
  vec_coeff->SetTime(t);
  gf->ProjectCoefficient(*vec_coeff);
}

CoupledCoefficient::CoupledCoefficient(
    const hephaestus::InputParameters &params)
    : AuxKernel(params),
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

} // namespace hephaestus

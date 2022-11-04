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

#include "kernels.hpp"

namespace hephaestus {

Kernel::Kernel(const hephaestus::InputParameters &params)
    : variable_name(params.GetParam<std::string>("VariableName")) {}

WeakCurlCurlKernel::WeakCurlCurlKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void WeakCurlCurlKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  u_ = variables.Get(variable_name);
  coef = domain_properties.scalar_property_map[coef_name];

  curlCurl = new mfem::ParBilinearForm(u_->ParFESpace());
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
  curlCurl->Assemble();
}

void WeakCurlCurlKernel ::ApplyKernel(mfem::ParLinearForm *lf) {
  curlCurl->AddMultTranspose(*u_, *lf, -1.0);
};

CurlCurlKernel::CurlCurlKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void CurlCurlKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map[coef_name];
}

void CurlCurlKernel ::ApplyKernel(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
};

VectorFEMassKernel::VectorFEMassKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void VectorFEMassKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map[coef_name];
}

void VectorFEMassKernel ::ApplyKernel(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*coef));
};

} // namespace hephaestus

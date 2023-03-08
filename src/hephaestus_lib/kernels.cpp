#include "kernels.hpp"

namespace hephaestus {

WeakCurlCurlKernel::WeakCurlCurlKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coupled_gf_name(params.GetParam<std::string>("CoupledVariableName")),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void WeakCurlCurlKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  u_ = variables.Get(coupled_gf_name);
  coef = domain_properties.scalar_property_map[coef_name];

  curlCurl = new mfem::ParBilinearForm(u_->ParFESpace());
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
  curlCurl->Assemble();
}

void WeakCurlCurlKernel ::Apply(mfem::ParLinearForm *lf) {
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

void CurlCurlKernel ::Apply(mfem::ParBilinearForm *blf) {
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

void VectorFEMassKernel ::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*coef));
};

DiffusionKernel::DiffusionKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void DiffusionKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map[coef_name];
}

void DiffusionKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::DiffusionIntegrator(*coef));
};

MixedVectorGradientKernel::MixedVectorGradientKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void MixedVectorGradientKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map[coef_name];
}

void MixedVectorGradientKernel::Apply(mfem::ParMixedBilinearForm *mblf) {
  mblf->AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(*coef));
};

VectorFEWeakDivergenceKernel::VectorFEWeakDivergenceKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void VectorFEWeakDivergenceKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map[coef_name];
}

void VectorFEWeakDivergenceKernel::Apply(mfem::ParMixedBilinearForm *mblf) {
  mblf->AddDomainIntegrator(new mfem::VectorFEWeakDivergenceIntegrator(*coef));
};

} // namespace hephaestus

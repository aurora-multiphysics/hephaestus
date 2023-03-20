#include "kernels.hpp"

namespace hephaestus {

WeakCurlKernel::WeakCurlKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      hcurl_gf_name(params.GetParam<std::string>("HCurlVarName")),
      hdiv_gf_name(params.GetParam<std::string>("HDivVarName")),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

WeakCurlKernel::~WeakCurlKernel() {
  if (weakCurl != NULL) {
    delete weakCurl;
  }
}

void WeakCurlKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  u_ = variables.Get(hcurl_gf_name);
  v_ = variables.Get(hdiv_gf_name);

  coef = domain_properties.scalar_property_map[coef_name];

  weakCurl = new mfem::ParMixedBilinearForm(u_->ParFESpace(), v_->ParFESpace());
  weakCurl->AddDomainIntegrator(new mfem::VectorFECurlIntegrator(*coef));
}

void WeakCurlKernel ::Apply(mfem::ParLinearForm *lf) {
  weakCurl->Update();
  weakCurl->Assemble();
  weakCurl->AddMultTranspose(*v_, *lf);
};

WeakCurlCurlKernel::WeakCurlCurlKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coupled_gf_name(params.GetParam<std::string>("CoupledVariableName")),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

WeakCurlCurlKernel::~WeakCurlCurlKernel() {
  if (curlCurl != NULL) {
    delete curlCurl;
  }
}

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

#include "weak_curl_curl_kernel.hpp"

namespace hephaestus {

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

void WeakCurlCurlKernel::Apply(mfem::ParLinearForm *lf) {
  curlCurl->AddMultTranspose(*u_, *lf, -1.0);
}

} // namespace hephaestus

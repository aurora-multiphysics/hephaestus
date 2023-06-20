#include "curl_curl_kernel.hpp"

namespace hephaestus {

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

void CurlCurlKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
}

} // namespace hephaestus

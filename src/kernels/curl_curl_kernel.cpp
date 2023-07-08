#include "curl_curl_kernel.hpp"

namespace hephaestus {

CurlCurlKernel::CurlCurlKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void CurlCurlKernel::Init(hephaestus::GridFunctions &gridfunctions,
                          const hephaestus::FESpaces &fespaces,
                          hephaestus::BCMap &bc_map,
                          hephaestus::Coefficients &coefficients) {

  coef = coefficients.scalars.Get(coef_name);
}

void CurlCurlKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
}

} // namespace hephaestus

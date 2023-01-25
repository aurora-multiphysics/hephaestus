#include "kernels.hpp"

namespace hephaestus {

WeakCurlCurlKernel::WeakCurlCurlKernel(
    const hephaestus::InputParameters &params)
    : gf_name(params.GetParam<std::string>("VariableName")),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void WeakCurlCurlKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  u_ = variables.Get(gf_name);
  coef = domain_properties.scalar_property_map[coef_name];

  curlCurl = new mfem::ParBilinearForm(u_->ParFESpace());
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
  curlCurl->Assemble();
}

void WeakCurlCurlKernel ::ApplyKernel(mfem::ParLinearForm *lf) {
  curlCurl->AddMultTranspose(*u_, *lf, -1.0);
};

} // namespace hephaestus

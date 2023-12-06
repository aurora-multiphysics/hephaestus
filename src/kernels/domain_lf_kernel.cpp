#include "domain_lf_kernel.hpp"

namespace hephaestus {

DomainLFKernel::DomainLFKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void DomainLFKernel::Init(hephaestus::GridFunctions &gridfunctions,
                           const hephaestus::FESpaces &fespaces,
                           hephaestus::BCMap &bc_map,
                           hephaestus::Coefficients &coefficients) {

  coef = coefficients.scalars.Get(coef_name);
}

void DomainLFKernel::Apply(mfem::ParLinearForm *lf) {
  lf->AddDomainIntegrator(new mfem::DomainLFIntegrator(*coef));
};

} // namespace hephaestus

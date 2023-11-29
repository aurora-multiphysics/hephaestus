#include "domain_divergence_lf_kernel.hpp"

namespace hephaestus {

DomainDivergenceLFKernel::DomainDivergenceLFKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void DomainDivergenceLFKernel::Init(hephaestus::GridFunctions &gridfunctions,
                          const hephaestus::FESpaces &fespaces,
                          hephaestus::BCMap &bc_map,
                          hephaestus::Coefficients &coefficients) {

  coef = coefficients.scalars.Get(coef_name);
}

void DomainDivergenceLFKernel::Apply(mfem::ParLinearForm *lf) {
  lf->AddDomainIntegrator(new mfem::DomainLFH1DivIntegrator(*coef));
}

} // namespace hephaestus

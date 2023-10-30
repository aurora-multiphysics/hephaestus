#include "mixed_weak_divergence_kernel.hpp"

namespace hephaestus {

MixedWeakDivergenceKernel::MixedWeakDivergenceKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void MixedWeakDivergenceKernel::Init(hephaestus::GridFunctions &gridfunctions,
                                     const hephaestus::FESpaces &fespaces,
                                     hephaestus::BCMap &bc_map,
                                     hephaestus::Coefficients &coefficients) {

  coef = coefficients.scalars.Get(coef_name);
}

void MixedWeakDivergenceKernel::Apply(mfem::ParMixedBilinearForm *mblf) {
  mblf->AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(*coef));
}

} // namespace hephaestus

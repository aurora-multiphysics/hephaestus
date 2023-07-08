#include "mixed_vector_gradient_kernel.hpp"

namespace hephaestus {

MixedVectorGradientKernel::MixedVectorGradientKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void MixedVectorGradientKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map, hephaestus::Coefficients &coefficients) {

  coef = coefficients.scalars.Get(coef_name);
}

void MixedVectorGradientKernel::Apply(mfem::ParMixedBilinearForm *mblf) {
  mblf->AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(*coef));
}

} // namespace hephaestus

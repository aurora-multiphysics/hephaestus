#include "diffusion_kernel.hpp"

namespace hephaestus {

DiffusionKernel::DiffusionKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void DiffusionKernel::Init(hephaestus::GridFunctions &gridfunctions,
                           const hephaestus::FESpaces &fespaces,
                           hephaestus::BCMap &bc_map,
                           hephaestus::Coefficients &coefficients) {

  coef = coefficients.scalars.Get(coef_name);
}

void DiffusionKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::DiffusionIntegrator(*coef));
};

} // namespace hephaestus

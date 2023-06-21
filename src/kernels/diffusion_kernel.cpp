#include "diffusion_kernel.hpp"

namespace hephaestus {

DiffusionKernel::DiffusionKernel(const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void DiffusionKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map.Get(coef_name);
}

void DiffusionKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::DiffusionIntegrator(*coef));
};

} // namespace hephaestus
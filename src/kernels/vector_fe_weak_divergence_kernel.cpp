#include "vector_fe_weak_divergence_kernel.hpp"

namespace hephaestus {

VectorFEWeakDivergenceKernel::VectorFEWeakDivergenceKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void VectorFEWeakDivergenceKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  coef = domain_properties.scalar_property_map[coef_name];
}

void VectorFEWeakDivergenceKernel::Apply(mfem::ParMixedBilinearForm *mblf) {
  mblf->AddDomainIntegrator(new mfem::VectorFEWeakDivergenceIntegrator(*coef));
};

} // namespace hephaestus

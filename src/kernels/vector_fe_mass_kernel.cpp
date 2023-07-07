#include "vector_fe_mass_kernel.hpp"

namespace hephaestus {

VectorFEMassKernel::VectorFEMassKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      coef_name(params.GetParam<std::string>("CoefficientName")) {}

void VectorFEMassKernel::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map, hephaestus::Coefficients &domain_properties) {

  coef = domain_properties.scalar_property_map.Get(coef_name);
}

void VectorFEMassKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*coef));
};

} // namespace hephaestus

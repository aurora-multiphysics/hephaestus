#include "vector_fe_mass_kernel.hpp"

namespace hephaestus
{

VectorFEMassKernel::VectorFEMassKernel(const hephaestus::InputParameters & params)
  : Kernel(params), coef_name(params.GetParam<std::string>("CoefficientName"))
{
}

void
VectorFEMassKernel::Init(hephaestus::GridFunctions & gridfunctions,
                         const hephaestus::FESpaces & fespaces,
                         hephaestus::BCMap & bc_map,
                         hephaestus::Coefficients & coefficients)
{

  coef = coefficients.scalars.Get(coef_name);
}

void
VectorFEMassKernel::Apply(mfem::ParBilinearForm * blf)
{
  blf->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*coef));
};

} // namespace hephaestus

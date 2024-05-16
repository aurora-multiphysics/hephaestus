#include "linear_kernel.hpp"

namespace hephaestus
{

LinearKernel::LinearKernel(const hephaestus::InputParameters & params)
  : Kernel(params), _coef_name(params.GetParam<std::string>("CoefficientName"))
{
}

void
LinearKernel::Init(hephaestus::GridFunctions & gridfunctions,
                   const hephaestus::FESpaces & fespaces,
                   hephaestus::BCMap & bc_map,
                   hephaestus::Coefficients & coefficients)
{
  _coef = coefficients._scalars.Get(_coef_name);
}

void
LinearKernel::Apply(mfem::ParLinearForm * lf)
{
  lf->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef));
}

} // namespace hephaestus

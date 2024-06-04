#include "weak_curl_kernel.hpp"

namespace hephaestus
{

WeakCurlKernel::WeakCurlKernel(const hephaestus::InputParameters & params)
  : Kernel(params),
    _hcurl_gf_name(params.Get<std::string>("HCurlVarName")),
    _hdiv_gf_name(params.Get<std::string>("HDivVarName")),
    _coef_name(params.Get<std::string>("CoefficientName"))
{
}

void
WeakCurlKernel::Init(hephaestus::GridFunctions & gridfunctions,
                     const hephaestus::FESpaces & fespaces,
                     hephaestus::BCMap & bc_map,
                     hephaestus::Coefficients & coefficients)
{
  _u = gridfunctions.Get(_hcurl_gf_name);
  _v = gridfunctions.Get(_hdiv_gf_name);

  _coef = coefficients._scalars.Get(_coef_name);

  _weak_curl = std::make_unique<mfem::ParMixedBilinearForm>(_u->ParFESpace(), _v->ParFESpace());
  _weak_curl->AddDomainIntegrator(new mfem::VectorFECurlIntegrator(*_coef));
  _weak_curl->Assemble();
}

void
WeakCurlKernel::Update()
{
  // NB: we only have ownership of _weak_curl! _u, _v will be updated by the
  // gridfunctions. Ensure gridfunctions are updated first.
  _weak_curl->Update();
  _weak_curl->Assemble();
}

void
WeakCurlKernel::Apply(mfem::ParLinearForm * lf)
{
  _weak_curl->AddMultTranspose(*_v, *lf);
}

} // namespace hephaestus

#include "weak_curl_kernel.hpp"

namespace hephaestus
{

WeakCurlKernel::WeakCurlKernel(const hephaestus::InputParameters & params)
  : Kernel(params),
    hcurl_gf_name(params.GetParam<std::string>("HCurlVarName")),
    hdiv_gf_name(params.GetParam<std::string>("HDivVarName")),
    coef_name(params.GetParam<std::string>("CoefficientName"))
{
}

WeakCurlKernel::~WeakCurlKernel()
{
  if (weakCurl != NULL)
  {
    delete weakCurl;
  }
}

void
WeakCurlKernel::Init(hephaestus::GridFunctions & gridfunctions,
                     const hephaestus::FESpaces & fespaces,
                     hephaestus::BCMap & bc_map,
                     hephaestus::Coefficients & coefficients)
{

  u_ = gridfunctions.Get(hcurl_gf_name);
  v_ = gridfunctions.Get(hdiv_gf_name);

  coef = coefficients.scalars.Get(coef_name);

  weakCurl = new mfem::ParMixedBilinearForm(u_->ParFESpace(), v_->ParFESpace());
  weakCurl->AddDomainIntegrator(new mfem::VectorFECurlIntegrator(*coef));
}

void
WeakCurlKernel::Apply(mfem::ParLinearForm * lf)
{
  weakCurl->Update();
  weakCurl->Assemble();
  weakCurl->AddMultTranspose(*v_, *lf);
}

} // namespace hephaestus

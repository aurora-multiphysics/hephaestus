#include "weak_curl_curl_kernel.hpp"

namespace hephaestus
{

WeakCurlCurlKernel::WeakCurlCurlKernel(const hephaestus::InputParameters & params)
  : Kernel(params),
    coupled_gf_name(params.GetParam<std::string>("CoupledVariableName")),
    coef_name(params.GetParam<std::string>("CoefficientName"))
{
}

WeakCurlCurlKernel::~WeakCurlCurlKernel()
{
  if (curlCurl != NULL)
  {
    delete curlCurl;
  }
}

void
WeakCurlCurlKernel::Init(hephaestus::GridFunctions & gridfunctions,
                         const hephaestus::FESpaces & fespaces,
                         hephaestus::BCMap & bc_map,
                         hephaestus::Coefficients & coefficients)
{

  u_ = gridfunctions.Get(coupled_gf_name);
  coef = coefficients.scalars.Get(coef_name);

  curlCurl = new mfem::ParBilinearForm(u_->ParFESpace());
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*coef));
  curlCurl->Assemble();
}

void
WeakCurlCurlKernel::Apply(mfem::ParLinearForm * lf)
{
  curlCurl->AddMultTranspose(*u_, *lf, -1.0);
}

} // namespace hephaestus

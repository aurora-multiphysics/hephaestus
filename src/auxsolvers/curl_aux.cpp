#include "curl_aux.hpp"

#include <utility>

namespace hephaestus
{

CurlAuxSolver::CurlAuxSolver(std::string input_gf_name, std::string curl_gf_name)
  : _input_gf_name(std::move(input_gf_name)), _curl_gf_name(std::move(curl_gf_name))
{
}

void
CurlAuxSolver::Init(const hephaestus::GridFunctions & gridfunctions,
                    hephaestus::Coefficients & coefficients)
{
  _u = gridfunctions.GetShared(_input_gf_name);
  if (_u == nullptr)
  {
    MFEM_ABORT("GridFunction " << _input_gf_name << " not found when initializing CurlAuxSolver");
  }
  _curl_u = gridfunctions.GetShared(_curl_gf_name);
  if (_curl_u == nullptr)
  {
    MFEM_ABORT("GridFunction " << _curl_gf_name << " not found when initializing CurlAuxSolver");
  }
  _curl =
      std::make_unique<mfem::ParDiscreteLinearOperator>(_u->ParFESpace(), _curl_u->ParFESpace());
  _curl->AddDomainInterpolator(new mfem::CurlInterpolator());
  _curl->Assemble();
}

void
CurlAuxSolver::Solve(double t)
{
  _curl->Mult(*_u, *_curl_u);
}

} // namespace hephaestus

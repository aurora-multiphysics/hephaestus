#include "curl_aux.hpp"

namespace hephaestus {

CurlAuxSolver::CurlAuxSolver(const std::string &input_gf_name,
                             const std::string &curl_gf_name)
    : AuxSolver(), _input_gf_name(input_gf_name), _curl_gf_name(curl_gf_name) {}

void CurlAuxSolver::Init(const hephaestus::GridFunctions &gridfunctions,
                         hephaestus::Coefficients &coefficients) {
  u_ = gridfunctions.Get(_input_gf_name);
  if (u_ == NULL) {
    MFEM_ABORT("GridFunction " << _input_gf_name
                               << " not found when initializing CurlAuxSolver");
  }
  curl_u_ = gridfunctions.Get(_curl_gf_name);
  if (curl_u_ == NULL) {
    MFEM_ABORT("GridFunction " << _curl_gf_name
                               << " not found when initializing CurlAuxSolver");
  }
  curl = std::make_unique<mfem::ParDiscreteLinearOperator>(
      u_->ParFESpace(), curl_u_->ParFESpace());
  curl->AddDomainInterpolator(new mfem::CurlInterpolator());
  curl->Assemble();
}

void CurlAuxSolver::Solve(double t) { curl->Mult(*u_, *curl_u_); }

} // namespace hephaestus

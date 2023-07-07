#include "curl_aux.hpp"

namespace hephaestus {

CurlAuxSolver::CurlAuxSolver(const hephaestus::InputParameters &params)
    : AuxSolver(), var_name(params.GetParam<std::string>("VariableName")),
      curl_var_name(params.GetParam<std::string>("CurlVariableName")) {}

void CurlAuxSolver::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::Coefficients &domain_properties) {
  u_ = variables.Get(var_name);
  curl_u_ = variables.Get(curl_var_name);

  curl = new mfem::ParDiscreteLinearOperator(u_->ParFESpace(),
                                             curl_u_->ParFESpace());
  curl->AddDomainInterpolator(new mfem::CurlInterpolator());
  curl->Assemble();
}

void CurlAuxSolver::Solve(double t) { curl->Mult(*u_, *curl_u_); }

} // namespace hephaestus

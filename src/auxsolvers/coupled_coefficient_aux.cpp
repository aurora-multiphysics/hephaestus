#include "coupled_coefficient_aux.hpp"

namespace hephaestus {

CoupledCoefficient::CoupledCoefficient(
    const hephaestus::InputParameters &params)
    : AuxSolver(),
      coupled_var_name(params.GetParam<std::string>("CoupledVariableName")) {}

void CoupledCoefficient::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::Coefficients &coefficients) {
  if (variables.Has(coupled_var_name)) {
    gf = variables.Get(coupled_var_name);
  } else {
    const std::string error_message = coupled_var_name +
                                      " not found in variables when "
                                      "creating CoupledCoefficient\n";
    mfem::mfem_error(error_message.c_str());
  }
}

double CoupledCoefficient::Eval(mfem::ElementTransformation &T,
                                const mfem::IntegrationPoint &ip) {
  return gf->GetValue(T, ip);
}

} // namespace hephaestus

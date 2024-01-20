#include "coupled_coefficient_aux.hpp"

namespace hephaestus
{

CoupledCoefficient::CoupledCoefficient(const hephaestus::InputParameters & params)
  : coupled_var_name(params.GetParam<std::string>("CoupledVariableName"))
{
}

void
CoupledCoefficient::Init(const hephaestus::GridFunctions & gridfunctions,
                         hephaestus::Coefficients & coefficients)
{
  if (gridfunctions.Has(coupled_var_name))
  {
    gf = gridfunctions.Get(coupled_var_name);
  }
  else
  {
    const std::string error_message = coupled_var_name + " not found in gridfunctions when "
                                                         "creating CoupledCoefficient\n";
    mfem::mfem_error(error_message.c_str());
  }
}

double
CoupledCoefficient::Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip)
{
  return gf->GetValue(T, ip);
}

} // namespace hephaestus

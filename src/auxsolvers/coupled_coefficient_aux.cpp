#include "coupled_coefficient_aux.hpp"

namespace hephaestus
{

CoupledCoefficient::CoupledCoefficient(const hephaestus::InputParameters & params)
  : _coupled_var_name(params.GetParam<std::string>("CoupledVariableName"))
{
}

void
CoupledCoefficient::Init(const hephaestus::GridFunctions & gridfunctions,
                         hephaestus::Coefficients & coefficients)
{
  if (gridfunctions.Has(_coupled_var_name))
  {
    _gf = gridfunctions.GetShared(_coupled_var_name);
  }
  else
  {
    const std::string error_message = _coupled_var_name + " not found in gridfunctions when "
                                                          "creating CoupledCoefficient\n";
    mfem::mfem_error(error_message.c_str());
  }
}

double
CoupledCoefficient::Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip)
{
  return _gf->GetValue(T, ip);
}

} // namespace hephaestus

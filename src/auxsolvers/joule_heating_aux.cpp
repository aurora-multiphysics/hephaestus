#include "joule_heating_aux.hpp"

namespace hephaestus {

double JouleHeatingCoefficient::Eval(mfem::ElementTransformation &T,
                                     const mfem::IntegrationPoint &ip) {
  double thisSigma;
  thisSigma = sigma->Eval(T, ip);

  mfem::Vector E_re;
  mfem::Vector E_im;

  E_gf_re->GetVectorValue(T, ip, E_re);
  if (E_gf_im == NULL) {
    return thisSigma * (E_re * E_re);
  } else {
    E_gf_im->GetVectorValue(T, ip, E_im);
    return 0.5 * thisSigma * (E_re * E_re + E_im * E_im);
  }
}

// VariableName: name of joule heating GridFunction to store solution
// CoefficientName = name of JouleHeatingCoefficient to save
// ElectricFieldName: name of gridfunction of the electric field
// ConductivityCoefName: name of conductivity coefficient
JouleHeatingAuxSolver::JouleHeatingAuxSolver(
    const hephaestus::InputParameters &params)
    : CoefficientAuxSolver(params),
      electric_field_name(params.GetParam<std::string>("ElectricFieldName")),
      conductivity_coef_name(
          params.GetParam<std::string>("ConductivityCoefName")),
      complex_field(params.GetOptionalParam<bool>("ComplexField", false)),
      sigma(nullptr), E_gf_re(nullptr), E_gf_im(nullptr) {}

void JouleHeatingAuxSolver::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {

  sigma = domain_properties.scalar_property_map.Get(conductivity_coef_name);
  if (sigma == NULL) {
    MFEM_ABORT("Conductivity coefficient not found for Joule heating");
  }
  if (complex_field) {
    E_gf_re = variables.Get(electric_field_name + "_real");
    E_gf_im = variables.Get(electric_field_name + "_imag");
    if (E_gf_re == NULL) {
      MFEM_ABORT("Electric field real component not found for Joule heating");
    }
    if (E_gf_re == NULL) {
      MFEM_ABORT(
          "Electric field imaginary component not found for Joule heating");
    }
  } else {
    E_gf_re = variables.Get(electric_field_name);
    if (E_gf_re == NULL) {
      MFEM_ABORT("Electric field not found for Joule heating");
    }
  }

  domain_properties.scalar_property_map.Register(
      coef_name, new JouleHeatingCoefficient(sigma, E_gf_re, E_gf_im), true);

  CoefficientAuxSolver::Init(variables, domain_properties);
}

} // namespace hephaestus

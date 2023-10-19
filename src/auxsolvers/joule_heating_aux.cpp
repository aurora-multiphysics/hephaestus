#include "joule_heating_aux.hpp"

namespace hephaestus {

double JouleHeatingDensityCoefficient::Eval(mfem::ElementTransformation &T,
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
// CoefficientName = name of LorentzForceCoefficient to save
// ElectricFieldName: name of gridfunction of the electric field
// ConductivityCoefName: name of conductivity coefficient
JouleHeatingDensityAux::JouleHeatingDensityAux(
    const std::string &joule_heating_density_gf_name,
    const std::string &joule_heating_density_coef_name,
    const std::string &electric_field_gf_name,
    const std::string &electric_conductivity_coef_name,
    const bool complex_field)
    : CoefficientAux(joule_heating_density_gf_name,
                     joule_heating_density_coef_name),
      _electric_field_gf_name(electric_field_gf_name),
      _electric_conductivity_coef_name(electric_conductivity_coef_name),
      _complex_field(complex_field), sigma(nullptr), E_gf_re(nullptr),
      E_gf_im(nullptr) {}

void JouleHeatingDensityAux::Init(
    const hephaestus::GridFunctions &gridfunctions,
    hephaestus::Coefficients &coefficients) {

  sigma = coefficients.scalars.Get(_electric_conductivity_coef_name);
  if (sigma == NULL) {
    MFEM_ABORT("Conductivity coefficient not found for Joule heating");
  }
  if (_complex_field) {
    E_gf_re = gridfunctions.Get(_electric_field_gf_name + "_real");
    E_gf_im = gridfunctions.Get(_electric_field_gf_name + "_imag");
    if (E_gf_re == NULL) {
      MFEM_ABORT("Electric field real component not found for Joule heating");
    }
    if (E_gf_im == NULL) {
      MFEM_ABORT(
          "Electric field imaginary component not found for Joule heating");
    }
  } else {
    E_gf_re = gridfunctions.Get(_electric_field_gf_name);
    if (E_gf_re == NULL) {
      MFEM_ABORT("Electric field not found for Joule heating");
    }
  }

  coefficients.scalars.Register(
      _coef_name, new JouleHeatingDensityCoefficient(sigma, E_gf_re, E_gf_im),
      true);

  CoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus

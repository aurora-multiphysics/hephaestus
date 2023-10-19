#pragma once
#include "auxsolvers.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus {

// The JouleHeatingDensityCoefficient object will contain a reference to the
// electric field grid function, and the conductivity sigma, and returns sigma E
// dot E at a point.
class JouleHeatingDensityCoefficient : public mfem::Coefficient {
private:
  mfem::ParGridFunction *E_gf_re;
  mfem::ParGridFunction *E_gf_im;
  mfem::Coefficient *sigma;

public:
  JouleHeatingDensityCoefficient(mfem::Coefficient *sigma_,
                                 mfem::ParGridFunction *E_gf_re_,
                                 mfem::ParGridFunction *E_gf_im_ = NULL)
      : E_gf_re(E_gf_re_), E_gf_im(E_gf_im_), sigma(sigma_) {}
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);
  virtual ~JouleHeatingDensityCoefficient() {}
};

// Project a stored scalar Coefficient onto a (scalar) GridFunction
class JouleHeatingDensityAux : public CoefficientAux {
private:
  mfem::Coefficient *sigma;
  mfem::ParGridFunction *E_gf_re;
  mfem::ParGridFunction *E_gf_im;

  std::string _electric_field_gf_name;
  std::string _electric_conductivity_coef_name;
  bool _complex_field;

public:
  JouleHeatingDensityAux(const std::string &joule_heating_density_gf_name,
                         const std::string &joule_heating_density_coef_name,
                         const std::string &electric_field_gf_name,
                         const std::string &electric_conductivity_coef_name,
                         const bool complex_field = false);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;
};

} // namespace hephaestus

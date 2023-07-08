#pragma once
#include "auxsolvers.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus {

// The JouleHeatingCoefficient object will contain a reference to the electric
// field grid function, and the conductivity sigma, and returns sigma E dot E at
// a point.
class JouleHeatingCoefficient : public mfem::Coefficient {
private:
  mfem::ParGridFunction *E_gf_re;
  mfem::ParGridFunction *E_gf_im;
  mfem::Coefficient *sigma;

public:
  JouleHeatingCoefficient(mfem::Coefficient *sigma_,
                          mfem::ParGridFunction *E_gf_re_,
                          mfem::ParGridFunction *E_gf_im_ = NULL)
      : E_gf_re(E_gf_re_), E_gf_im(E_gf_im_), sigma(sigma_) {}
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);
  virtual ~JouleHeatingCoefficient() {}
};

// Project a stored scalar Coefficient onto a (scalar) GridFunction
class JouleHeatingAuxSolver : public CoefficientAuxSolver {
private:
  mfem::Coefficient *sigma;
  mfem::ParGridFunction *E_gf_re;
  mfem::ParGridFunction *E_gf_im;

  std::string electric_field_name;
  std::string conductivity_coef_name;
  bool complex_field;

public:
  JouleHeatingAuxSolver(const hephaestus::InputParameters &params);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;
};

} // namespace hephaestus

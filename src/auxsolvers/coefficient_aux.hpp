#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Project a stored scalar Coefficient onto a (scalar) GridFunction
class CoefficientAux : public AuxSolver {
public:
  CoefficientAux(const std::string &gf_name, const std::string &coef_name);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;

  void Solve(double t = 0.0) override;

protected:
  const std::string _gf_name;   // name of the variable
  const std::string _coef_name; // name of the coefficient

  mfem::ParGridFunction *gf;
  mfem::Coefficient *coef;
};

} // namespace hephaestus

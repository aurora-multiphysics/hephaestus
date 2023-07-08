#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Project a stored scalar Coefficient onto a (scalar) GridFunction
class CoefficientAuxSolver : public AuxSolver {
public:
  CoefficientAuxSolver(const hephaestus::InputParameters &params);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;

  void Solve(double t = 0.0) override;

  std::string var_name;  // name of the variable
  std::string coef_name; // name of the coefficient

  mfem::ParGridFunction *gf;
  mfem::Coefficient *coeff;
};

} // namespace hephaestus

#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus
{

// Project a stored scalar Coefficient onto a (scalar) GridFunction
class CoefficientAux : public AuxSolver
{
public:
  CoefficientAux(std::string gf_name, std::string coef_name);

  ~CoefficientAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  void Solve(double t = 0.0) override;

protected:
  const std::string _gf_name;   // name of the variable
  const std::string _coef_name; // name of the coefficient

  mfem::ParGridFunction * _gf{nullptr};
  mfem::Coefficient * _coef{nullptr};
};

} // namespace hephaestus

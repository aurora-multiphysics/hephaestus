#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Project a stored vector Coefficient onto a (vector) GridFunction
class VectorCoefficientAux : public AuxSolver {
public:
  VectorCoefficientAux(const std::string &gf_name,
                       const std::string &vec_coef_name);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;

  void Solve(double t = 0.0) override;

protected:
  const std::string _gf_name;       // name of the variable
  const std::string _vec_coef_name; // name of the vector coefficient

  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coef;
};

} // namespace hephaestus

#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Project a stored vector Coefficient onto a (vector) GridFunction
class VectorCoefficientAuxSolver : public AuxSolver {
public:
  VectorCoefficientAuxSolver(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  void Solve(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string vec_coef_name; // name of the vector coefficient

  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coeff;
};

} // namespace hephaestus

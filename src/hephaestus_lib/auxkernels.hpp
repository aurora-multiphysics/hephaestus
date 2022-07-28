#pragma once
#include "materials.hpp"
#include "mfem.hpp"
#include "variables.hpp"

namespace hephaestus {

class VectorCoefficientAuxKernel {
public:
  VectorCoefficientAuxKernel(const std::string &var_name_,
                             const std::string &vec_coef_name_);

  void Init(const hephaestus::VariableMap &variables,
            hephaestus::DomainProperties &domain_properties);

  void Solve(double t);

  std::string var_name;      // name of the variable
  std::string vec_coef_name; // name of the vector coefficient

  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coeff;
};

} // namespace hephaestus

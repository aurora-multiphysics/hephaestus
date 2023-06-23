#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Calculate the curl of a gridfunction.
class CurlAuxSolver : public AuxSolver {
public:
  CurlAuxSolver(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  void Solve(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string curl_var_name; // Variable in which to store curl

  mfem::ParGridFunction *u_, *curl_u_;
  mfem::ParDiscreteLinearOperator *curl;
};

} // namespace hephaestus

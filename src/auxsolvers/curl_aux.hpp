#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus
{

// Calculate the curl of a gridfunction.
class CurlAuxSolver : public AuxSolver
{
public:
  CurlAuxSolver(const std::string & input_gf_name, const std::string & curl_gf_name);

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  void Solve(double t = 0.0) override;

private:
  const std::string _input_gf_name; // name of the variable
  const std::string _curl_gf_name;  // Variable in which to store curl

  mfem::ParGridFunction *u_, *curl_u_;
  mfem::ParDiscreteLinearOperator * curl;
};

} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class StaticsFormulation : public SteadyStateEMFormulation {
public:
  StaticsFormulation(const std::string &alpha_coef_name,
                     const std::string &h_curl_var_name);

  virtual void ConstructJacobianPreconditioner() override;

  virtual void ConstructJacobianSolver() override;

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _h_curl_var_name;
};

class StaticsOperator : public ProblemOperator {
public:
  StaticsOperator(mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
                  hephaestus::GridFunctions &gridfunctions,
                  hephaestus::BCMap &bc_map,
                  hephaestus::Coefficients &coefficients,
                  hephaestus::Sources &sources, mfem::Solver &jacobian_solver,
                  const std::string &h_curl_var_name,
                  const std::string &stiffness_coef_name);

  ~StaticsOperator(){};

  virtual void SetGridFunctions() override;
  virtual void Init(mfem::Vector &X) override;
  virtual void Solve(mfem::Vector &X) override;

private:
  const std::string _h_curl_var_name, _stiffness_coef_name;

  mfem::Coefficient *stiffCoef_; // Stiffness Material Coefficient
};

} // namespace hephaestus

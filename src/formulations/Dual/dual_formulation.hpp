#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class DualFormulation : public TimeDomainEMFormulation {
public:
  DualFormulation(const std::string &alpha_coef_name,
                  const std::string &beta_coef_name,
                  const std::string &h_curl_var_name,
                  const std::string &h_div_var_name);

  virtual void ConstructEquationSystem() override;

  virtual void ConstructJacobianPreconditioner() override;

  virtual void ConstructJacobianSolver() override;

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _h_curl_var_name;
  const std::string _h_div_var_name;
};

class WeakCurlEquationSystem : public TimeDependentEquationSystem {
public:
  WeakCurlEquationSystem(const hephaestus::InputParameters &params);

  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void addKernels() override;

  std::string _h_curl_var_name, _h_div_var_name, _alpha_coef_name,
      _beta_coef_name, _dtalpha_coef_name;
};

class DualOperator : public TimeDomainProblemOperator {
public:
  DualOperator(hephaestus::Problem &problem)
      : TimeDomainProblemOperator(problem){};

  ~DualOperator(){};

  void Init(mfem::Vector &X) override;

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;
  virtual void SetGridFunctions() override;

  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParFiniteElementSpace *HDivFESpace_;

  std::string _h_curl_var_name, _h_div_var_name;

  mfem::ParGridFunction *u_;  // HCurl vector field
  mfem::ParGridFunction *dv_; // HDiv vector field

protected:
  mfem::ParDiscreteLinearOperator *curl;
};
} // namespace hephaestus

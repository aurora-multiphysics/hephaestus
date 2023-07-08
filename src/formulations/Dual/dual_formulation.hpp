#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class DualFormulation : public TimeDomainFormulation {
public:
  DualFormulation();

  virtual void ConstructEquationSystem() override;

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;

protected:
  std::string h_curl_var_name, h_div_var_name, alpha_coef_name, beta_coef_name;
};

class WeakCurlEquationSystem : public TimeDependentEquationSystem {
public:
  WeakCurlEquationSystem(const hephaestus::InputParameters &params);

  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void addKernels() override;

  std::string h_curl_var_name, h_div_var_name, alpha_coef_name, beta_coef_name,
      dtalpha_coef_name;
};

class DualOperator : public TimeDomainEquationSystemOperator {
public:
  DualOperator(mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
               hephaestus::GridFunctions &gridfunctions,
               hephaestus::BCMap &bc_map,
               hephaestus::Coefficients &coefficients,
               hephaestus::Sources &sources,
               hephaestus::InputParameters &solver_options);

  ~DualOperator(){};

  void Init(mfem::Vector &X) override;

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;
  virtual void Setgridfunctions() override;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParFiniteElementSpace *HDivFESpace_;

  std::string h_curl_var_name, h_div_var_name;

  mfem::ParGridFunction *u_;  // HCurl vector field
  mfem::ParGridFunction *dv_; // HDiv vector field

protected:
  mfem::ParDiscreteLinearOperator *curl;
};
} // namespace hephaestus

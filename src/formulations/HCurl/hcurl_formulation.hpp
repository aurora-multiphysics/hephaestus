#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus
{

class HCurlFormulation : public TimeDomainEMFormulation
{
public:
  HCurlFormulation(const std::string & alpha_coef_name,
                   const std::string & beta_coef_name,
                   const std::string & h_curl_var_name);

  ~HCurlFormulation() override {}

  virtual void ConstructEquationSystem() override;

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _h_curl_var_name;
};

class CurlCurlEquationSystem : public TimeDependentEquationSystem
{
public:
  CurlCurlEquationSystem(const hephaestus::InputParameters & params);

  virtual void Init(hephaestus::GridFunctions & gridfunctions,
                    const hephaestus::FESpaces & fespaces,
                    hephaestus::BCMap & bc_map,
                    hephaestus::Coefficients & coefficients) override;
  virtual void addKernels() override;

  std::string h_curl_var_name, alpha_coef_name, beta_coef_name, dtalpha_coef_name;
};

class HCurlOperator : public TimeDomainEquationSystemOperator
{
public:
  HCurlOperator(mfem::ParMesh & pmesh,
                hephaestus::FESpaces & fespaces,
                hephaestus::GridFunctions & gridfunctions,
                hephaestus::BCMap & bc_map,
                hephaestus::Coefficients & coefficients,
                hephaestus::Sources & sources,
                hephaestus::InputParameters & solver_options);

  ~HCurlOperator() override {}

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
};

} // namespace hephaestus

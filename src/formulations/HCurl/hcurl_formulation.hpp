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
  HCurlFormulation(std::string alpha_coef_name,
                   std::string beta_coef_name,
                   std::string h_curl_var_name);

  ~HCurlFormulation() override = default;

  void ConstructEquationSystem() override;

  void ConstructOperator() override;

  void RegisterGridFunctions() override;

  void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _h_curl_var_name;
};

class CurlCurlEquationSystem : public TimeDependentEquationSystem
{
public:
  CurlCurlEquationSystem(const hephaestus::InputParameters & params);

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void AddKernels() override;

  std::string _h_curl_var_name, _alpha_coef_name, _beta_coef_name, _dtalpha_coef_name;
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
                hephaestus::InputParameters & solver_options,
                hephaestus::ProblemSolvers & solvers);

  ~HCurlOperator() override = default;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;

private:
  ProblemSolvers * _solvers;
};

} // namespace hephaestus

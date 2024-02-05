#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus
{

class DualFormulation : public TimeDomainEMFormulation
{
public:
  DualFormulation(std::string alpha_coef_name,
                  std::string beta_coef_name,
                  std::string h_curl_var_name,
                  std::string h_div_var_name);

  ~DualFormulation() override = default;

  void ConstructEquationSystem() override;

  void ConstructOperator() override;

  void RegisterGridFunctions() override;

  void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _h_curl_var_name;
  const std::string _h_div_var_name;
};

class WeakCurlEquationSystem : public TimeDependentEquationSystem
{
public:
  WeakCurlEquationSystem(const hephaestus::InputParameters & params);
  ~WeakCurlEquationSystem() override = default;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void AddKernels() override;

  std::string _h_curl_var_name, _h_div_var_name, _alpha_coef_name, _beta_coef_name,
      _dtalpha_coef_name;
};

class DualOperator : public TimeDomainEquationSystemOperator
{
public:
  DualOperator(mfem::ParMesh & pmesh,
               hephaestus::FESpaces & fespaces,
               hephaestus::GridFunctions & gridfunctions,
               hephaestus::BCMap & bc_map,
               hephaestus::Coefficients & coefficients,
               hephaestus::Sources & sources,
               hephaestus::InputParameters & solver_options,
               hephaestus::ProblemSolvers & solvers);

  ~DualOperator() override = default;

  void Init(mfem::Vector & X) override;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
  void SetGridFunctions() override;
  mfem::ParFiniteElementSpace * _h_curl_fe_space;
  mfem::ParFiniteElementSpace * _h_div_fe_space;

  std::string _h_curl_var_name, _h_div_var_name;

  mfem::ParGridFunction * _u;  // HCurl vector field
  mfem::ParGridFunction * _dv; // HDiv vector field

private:
  ProblemSolvers * _solvers;

protected:
  std::unique_ptr<mfem::ParDiscreteLinearOperator> _curl;
};
} // namespace hephaestus

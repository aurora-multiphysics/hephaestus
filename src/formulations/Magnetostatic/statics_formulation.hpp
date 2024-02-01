#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus
{

class StaticsFormulation : public SteadyStateEMFormulation
{
public:
  StaticsFormulation(std::string alpha_coef_name, std::string h_curl_var_name);

  ~StaticsFormulation() override = default;

  void ConstructOperator() override;

  void RegisterGridFunctions() override;

  void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _h_curl_var_name;
};

class StaticsOperator : public EquationSystemOperator
{
public:
  StaticsOperator(mfem::ParMesh & pmesh,
                  hephaestus::FESpaces & fespaces,
                  hephaestus::GridFunctions & gridfunctions,
                  hephaestus::BCMap & bc_map,
                  hephaestus::Coefficients & coefficients,
                  hephaestus::Sources & sources,
                  hephaestus::InputParameters & solver_options,
                  hephaestus::ProblemSolvers & solvers);

  ~StaticsOperator() override = default;

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;
  void Solve(mfem::Vector & X) override;

private:
  std::string _h_curl_var_name, _stiffness_coef_name;
  ProblemSolvers * _solvers;

  mfem::Coefficient * _stiff_coef; // Stiffness Material Coefficient
};

} // namespace hephaestus

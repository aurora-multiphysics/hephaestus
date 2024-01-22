#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus
{

class AVFormulation : public TimeDomainEMFormulation
{
public:
  AVFormulation(std::string alpha_coef_name,
                std::string inv_alpha_coef_name,
                std::string beta_coef_name,
                std::string vector_potential_name,
                std::string scalar_potential_name);

  ~AVFormulation() override = default;

  void ConstructEquationSystem() override;

  void ConstructOperator() override;

  void RegisterGridFunctions() override;

  void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _inv_alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _vector_potential_name;
  const std::string _scalar_potential_name;
};

class AVEquationSystem : public TimeDependentEquationSystem
{
public:
  AVEquationSystem(const hephaestus::InputParameters & params);

  ~AVEquationSystem() override = default;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void AddKernels() override;

  std::string a_name, v_name, coupled_variable_name, alpha_coef_name, beta_coef_name,
      dtalpha_coef_name, neg_beta_coef_name;
  mfem::ConstantCoefficient negCoef;
};

class AVOperator : public TimeDomainEquationSystemOperator
{
public:
  AVOperator(mfem::ParMesh & pmesh,
             hephaestus::FESpaces & fespaces,
             hephaestus::GridFunctions & gridfunctions,
             hephaestus::BCMap & bc_map,
             hephaestus::Coefficients & coefficients,
             hephaestus::Sources & sources,
             hephaestus::InputParameters & solver_options);

  ~AVOperator() override = default;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
};
} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class AVFormulation : public TimeDomainFormulation {
  std::string vector_potential_name, scalar_potential_name, alpha_coef_name,
      beta_coef_name;

public:
  AVFormulation();

  virtual void ConstructEquationSystem() override;

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;
};

class AVEquationSystem : public TimeDependentEquationSystem {
public:
  AVEquationSystem(const hephaestus::InputParameters &params);

  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients);
  virtual void addKernels() override;

  std::string a_name, v_name, coupled_variable_name, alpha_coef_name,
      beta_coef_name, dtalpha_coef_name, neg_beta_coef_name;
  mfem::ConstantCoefficient negCoef;
};

class AVOperator : public TimeDomainEquationSystemOperator {
public:
  AVOperator(mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
             hephaestus::GridFunctions &gridfunctions,
             hephaestus::BCMap &bc_map, hephaestus::Coefficients &coefficients,
             hephaestus::Sources &sources,
             hephaestus::InputParameters &solver_options);

  ~AVOperator(){};

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;
};
} // namespace hephaestus

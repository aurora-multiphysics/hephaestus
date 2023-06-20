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

  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties);
  virtual void addKernels() override;

  std::string a_name, v_name, coupled_variable_name, alpha_coef_name,
      beta_coef_name, dtalpha_coef_name, neg_beta_coef_name;
  mfem::ConstantCoefficient negCoef;
};

class AVOperator : public TimeDomainEquationSystemOperator {
public:
  AVOperator(mfem::ParMesh &pmesh,
             mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
             mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
             hephaestus::BCMap &bc_map,
             hephaestus::DomainProperties &domain_properties,
             hephaestus::Sources &sources,
             hephaestus::InputParameters &solver_options);

  ~AVOperator(){};

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;
};
} // namespace hephaestus

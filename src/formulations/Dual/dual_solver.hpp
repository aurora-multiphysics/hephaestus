#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class DualFormulation : public TransientFormulation {
public:
  DualFormulation();

  virtual hephaestus::TimeDependentEquationSystem *
  CreateEquationSystem() override;

  virtual hephaestus::TimeDomainEquationSystemOperator *
  CreateTimeDomainOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) override;

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) override;

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels) override{};

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;

protected:
  std::string h_curl_var_name, h_div_var_name, alpha_coef_name, beta_coef_name;
};

class DualOperator : public TimeDomainEquationSystemOperator {
public:
  DualOperator(mfem::ParMesh &pmesh,
               mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
               mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
               hephaestus::BCMap &bc_map,
               hephaestus::DomainProperties &domain_properties,
               hephaestus::Sources &sources,
               hephaestus::InputParameters &solver_options);

  ~DualOperator(){};

  void Init(mfem::Vector &X) override;

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;
  virtual void SetVariables() override;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParFiniteElementSpace *HDivFESpace_;

  std::string h_curl_var_name, h_div_var_name;

  mfem::ParGridFunction *u_;  // HCurl vector field
  mfem::ParGridFunction *dv_; // HDiv vector field

protected:
  mfem::ParDiscreteLinearOperator *curl;
};
} // namespace hephaestus

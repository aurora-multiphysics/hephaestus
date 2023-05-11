#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class AVFormulation : public TransientFormulation {
  std::string vector_potential_name, scalar_potential_name, alpha_coef_name,
      beta_coef_name;

public:
  AVFormulation();

  virtual std::unique_ptr<hephaestus::TimeDependentEquationSystem>
  CreateTimeDependentEquationSystem() const override;

  virtual std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator>
  CreateTimeDomainEquationSystemOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) const override;

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) override;

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels) override{};

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
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

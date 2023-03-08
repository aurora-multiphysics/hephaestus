#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HJDualFormulation : public hephaestus::DualFormulation {
public:
  HJDualFormulation();

  ~HJDualFormulation(){};

  virtual hephaestus::TimeDependentEquationSystem *
  CreateEquationSystem() override;

  virtual hephaestus::TimeDomainEquationSystemOperator *
  CreateTimeDomainOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) override;

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};

class HJDualOperator : public hephaestus::DualOperator {
public:
  HJDualOperator(mfem::ParMesh &pmesh, int order,
                 mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
                 mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                 hephaestus::BCMap &bc_map,
                 hephaestus::DomainProperties &domain_properties,
                 hephaestus::Sources &sources,
                 hephaestus::InputParameters &solver_options);

  ~HJDualOperator(){};

  virtual void SetMaterialCoefficients(
      hephaestus::DomainProperties &domain_properties) override;

  virtual void RegisterVariables() override;
};

} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EBDualSolver : public hephaestus::DualSolver {
  virtual void SetMaterialCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
  virtual void RegisterVariables() override;

public:
  EBDualSolver(mfem::ParMesh &pmesh, int order,
               mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
               hephaestus::BCMap &bc_map,
               hephaestus::DomainProperties &domain_properties);

  ~EBDualSolver(){};

  virtual void WriteConsoleSummary(double t, int it) override;
  double ElectricLosses() const;
};
} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HJDualSolver : public hephaestus::DualSolver {
  virtual void SetMaterialCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
  virtual void SetVariableNames() override;

public:
  HJDualSolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
               hephaestus::DomainProperties &domain_properties);

  ~HJDualSolver(){};
};
} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EFormSolver : public hephaestus::HCurlSolver {
  virtual void SetMaterialCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
  virtual void SetVariableNames() override;

public:
  EFormSolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
              hephaestus::DomainProperties &domain_properties);

  ~EFormSolver(){};
};
} // namespace hephaestus

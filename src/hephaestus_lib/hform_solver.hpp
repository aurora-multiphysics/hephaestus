#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HFormSolver : public hephaestus::HCurlSolver {
  virtual void SetMaterialCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
  virtual void SetVariableNames() override;

public:
  HFormSolver(mfem::ParMesh &pmesh, int order,
              mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
              mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
              hephaestus::BCMap &bc_map,
              hephaestus::DomainProperties &domain_properties);

  ~HFormSolver(){};
};
} // namespace hephaestus

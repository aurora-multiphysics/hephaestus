#pragma once
#include "../common/pfem_extras.hpp"
#include "eb_dual_solver.hpp"
#include "eform_solver.hpp"
#include "hform_solver.hpp"
#include "hj_dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class Factory {
public:
  static hephaestus::TransientFormulation *createTransientFormulation(
      std::string &formulation, mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties);

  static mfem::ParFiniteElementSpace *
  createParFESpace(const hephaestus::InputParameters params,
                   mfem::ParMesh &pmesh);
};

} // namespace hephaestus

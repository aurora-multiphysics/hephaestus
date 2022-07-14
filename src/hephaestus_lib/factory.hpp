#pragma once
#include "../common/pfem_extras.hpp"
#include "eb_dual_solver.hpp"
#include "hform_solver.hpp"
#include "hj_dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class FormulationFactory {
public:
  static hephaestus::TransientFormulation *
  createTransientFormulation(std::string &formulation, mfem::ParMesh &pmesh,
                             int order, hephaestus::BCMap &bc_map,
                             hephaestus::DomainProperties &domain_properties);
};

} // namespace hephaestus

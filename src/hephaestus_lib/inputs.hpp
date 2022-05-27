#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "boundary_conditions.hpp"
#include "executioner.hpp"
#include "materials.hpp"
#include "outputs.hpp"

namespace hephaestus {

class Inputs {
public:
  Inputs(){};
  Inputs(const mfem::Mesh &mesh_, const std::string &formulation_,
         const int order_, const BCMap &bc_map_,
         const DomainProperties &domain_properties_,
         const Executioner &executioner_, Outputs outputs_);

  mfem::Mesh mesh;
  std::string formulation;
  int order;
  BCMap bc_map;
  DomainProperties domain_properties;
  Executioner executioner;
  Outputs outputs;
};

} // namespace hephaestus
#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "boundary_conditions.hpp"
#include "executioner.hpp"
// #include "joule_solver.hpp"
#include "materials.hpp"

namespace hephaestus {

class Inputs {
public:
  Inputs(){};
  Inputs(const std::string &mesh_file_, const std::string &formulation_,
         const int order_, const BCMap &bc_map_,
         const DomainProperties &domain_properties_,
         const Executioner &executioner_);

  std::string mesh_file;
  std::string formulation;
  int order;
  BCMap bc_map;
  DomainProperties domain_properties;
  Executioner executioner;
};

} // namespace hephaestus
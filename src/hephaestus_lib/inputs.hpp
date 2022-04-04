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
  Inputs(const std::string &mesh_file, const std::string &formulation,
         const int order, const BCMap &boundary_conditions,
         const DomainProperties &domain_properties_,
         const Executioner &executioner_);

  std::string _mesh_file;
  std::string _formulation;
  int _order;
  BCMap bc_map;
  DomainProperties domain_properties;
  Executioner executioner;
};

} // namespace hephaestus
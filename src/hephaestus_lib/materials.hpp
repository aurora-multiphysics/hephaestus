#pragma once
#include <fstream>
#include <iostream>
#include <memory>

// #include "joule_solver.hpp"
#include "mesh_extras.hpp"

namespace hephaestus {

class Subdomain {
public:
  Subdomain(const std::string &subdomain_name, int subdomain_id);

  std::string name;
  int id;
  std::map<std::string, mfem::Coefficient *> property_map;
};

class DomainProperties {
public:
  DomainProperties();
  DomainProperties(std::vector<Subdomain> subdomains);

  mfem::PWCoefficient getGlobalScalarProperty(std::string property_name);

  std::map<std::string, mfem::Coefficient *> scalar_property_map;
  std::vector<Subdomain> _subdomains;
};

} // namespace hephaestus
#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

class Subdomain {
public:
  Subdomain(const std::string &name_, int id_);

  std::string name;
  int id;
  std::map<std::string, mfem::Coefficient *> property_map;
};

class DomainProperties {
public:
  DomainProperties();
  DomainProperties(std::vector<Subdomain> subdomains_);

  mfem::PWCoefficient getGlobalScalarProperty(std::string property_name_);

  std::map<std::string, mfem::Coefficient *> scalar_property_map;
  std::vector<Subdomain> subdomains;
};

} // namespace hephaestus
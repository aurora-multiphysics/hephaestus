#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

double prodFunc(double a, double b);
double fracFunc(double a, double b);

class Subdomain {
public:
  Subdomain(const std::string &name_, int id_);

  std::string name;
  int id;
  mfem::NamedFieldsMap<mfem::Coefficient> property_map;
};

class DomainProperties {
  double t; // Time at which time-dependent coefficients are evaluated
public:
  DomainProperties();
  DomainProperties(std::vector<Subdomain> subdomains_);
  void SetTime(double t);
  mfem::PWCoefficient getGlobalScalarProperty(std::string property_name_);

  mfem::NamedFieldsMap<mfem::Coefficient> scalar_property_map;
  mfem::NamedFieldsMap<mfem::VectorCoefficient> vector_property_map;
  std::vector<Subdomain> subdomains;
};

} // namespace hephaestus

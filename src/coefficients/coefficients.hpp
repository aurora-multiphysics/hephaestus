#pragma once
#include "mesh_extras.hpp"
#include "named_fields_map.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_set>

namespace hephaestus
{

double prodFunc(double a, double b);
double fracFunc(double a, double b);

class Subdomain
{
public:
  Subdomain(std::string name_, int id_);

  std::string name;
  int id;
  hephaestus::NamedFieldsMap<mfem::Coefficient> scalar_coefficients;
};

// Coefficients - stores all scalar and vector coefficients
//--SetTime
//--scalars
//--vectors

// Stores all coefficients defined over
class Coefficients
{
  double t; // Time at which time-dependent coefficients are evaluated
public:
  Coefficients();
  ~Coefficients() = default;

  Coefficients(std::vector<Subdomain> subdomains_);
  void SetTime(double t);
  void AddGlobalCoefficientsFromSubdomains();
  void registerDefaultCoefficients();

  hephaestus::NamedFieldsMap<mfem::Coefficient> scalars;
  hephaestus::NamedFieldsMap<mfem::VectorCoefficient> vectors;
  std::vector<Subdomain> subdomains;
};

} // namespace hephaestus

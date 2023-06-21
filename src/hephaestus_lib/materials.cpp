#include "materials.hpp"

namespace hephaestus {

double prodFunc(double a, double b) { return a * b; }
double fracFunc(double a, double b) { return a / b; }

Subdomain::Subdomain(const std::string &name_, int id_)
    : name(name_), id(id_) {}

DomainProperties::DomainProperties() {}

DomainProperties::DomainProperties(std::vector<Subdomain> subdomains_)
    : subdomains(subdomains_) {}

void DomainProperties::SetTime(double time) {
  for (auto const &[name, coeff_] : scalar_property_map) {
    coeff_->SetTime(time);
  }
  for (auto const &[name, vec_coeff_] : vector_property_map) {
    vec_coeff_->SetTime(time);
  }
  t = time;
}

mfem::PWCoefficient
DomainProperties::getGlobalScalarProperty(std::string property_name_) {

  mfem::Array<int> subdomain_ids;
  mfem::Array<mfem::Coefficient *> subdomain_coefs;

  for (std::size_t i = 0; i < subdomains.size(); i++) {
    subdomain_ids.Append(subdomains[i].id);
    subdomain_coefs.Append(subdomains[i].property_map.Get(property_name_));
  }
  mfem::PWCoefficient global_property_map(subdomain_ids, subdomain_coefs);
  return global_property_map;
}

} // namespace hephaestus

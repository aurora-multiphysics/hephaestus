#include "materials.hpp"

namespace hephaestus {

Subdomain::Subdomain(const std::string &name_, int id_)
    : name(name_), id(id_) {}

DomainProperties::DomainProperties() {}

DomainProperties::DomainProperties(std::vector<Subdomain> subdomains_)
    : subdomains(subdomains_) {}

mfem::PWCoefficient
DomainProperties::getGlobalScalarProperty(std::string property_name_) {

  mfem::Array<int> subdomain_ids;
  mfem::Array<mfem::Coefficient *> subdomain_coefs;

  for (std::size_t i = 0; i < subdomains.size(); i++) {
    subdomain_ids.Append(subdomains[i].id);
    subdomain_coefs.Append(subdomains[i].property_map[property_name_]);
  }
  mfem::PWCoefficient global_property_map(subdomain_ids, subdomain_coefs);
  return global_property_map;
}

} // namespace hephaestus
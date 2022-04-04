#include "materials.hpp"

namespace hephaestus {

Subdomain::Subdomain(const std::string &subdomain_name, int subdomain_id)
    : name(subdomain_name), id(subdomain_id) {}

DomainProperties::DomainProperties() {}

DomainProperties::DomainProperties(std::vector<Subdomain> subdomains)
    : _subdomains(subdomains) {}

mfem::PWCoefficient
DomainProperties::getGlobalScalarProperty(std::string property_name) {

  mfem::Array<int> subdomain_ids;
  mfem::Array<mfem::Coefficient *> subdomain_coefs;

  for (std::size_t i = 0; i < _subdomains.size(); i++) {
    subdomain_ids.Append(_subdomains[i].id);
    subdomain_coefs.Append(_subdomains[i].property_map[property_name]);
  }
  mfem::PWCoefficient global_property_map(subdomain_ids, subdomain_coefs);
  return global_property_map;
}

} // namespace hephaestus
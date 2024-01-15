#include "coefficients.hpp"

namespace hephaestus
{

double
prodFunc(double a, double b)
{
  return a * b;
}
double
fracFunc(double a, double b)
{
  return a / b;
}

Subdomain::Subdomain(const std::string & name_, int id_) : name(name_), id(id_) {}

Coefficients::Coefficients() { registerDefaultCoefficients(); }

Coefficients::Coefficients(std::vector<Subdomain> subdomains_) : subdomains(subdomains_)
{
  AddGlobalCoefficientsFromSubdomains();
  registerDefaultCoefficients();
}

void
Coefficients::registerDefaultCoefficients()
{
  scalars.Register("_one", new mfem::ConstantCoefficient(1.0), true);
}

void
Coefficients::SetTime(double time)
{
  for (auto const & [name, coeff_] : scalars)
  {
    coeff_->SetTime(time);
  }
  for (auto const & [name, vec_coeff_] : vectors)
  {
    vec_coeff_->SetTime(time);
  }
  t = time;
}

// merge subdomains?
void
Coefficients::AddGlobalCoefficientsFromSubdomains()
{

  // iterate over subdomains
  // check IDs span domain
  // accumulate list of property_name in unordered map

  // for each property_name
  // iterate over subdomains
  // accumulate coefs
  mfem::Array<int> subdomain_ids;
  std::unordered_set<std::string> scalar_property_names;

  for (std::size_t i = 0; i < subdomains.size(); i++)
  {
    subdomain_ids.Append(subdomains[i].id);
    // accumulate property names on subdomains, ignoring duplicates
    for (auto const & [name, coeff_] : subdomains[i].scalar_coefficients)
    {
      scalar_property_names.insert(name);
    }
  }
  // check if IDs span
  // iterate over properties stored on subdomains, and create global
  // coefficients
  for (auto & scalar_property_name : scalar_property_names)
  {
    mfem::Array<mfem::Coefficient *> subdomain_coefs;
    for (std::size_t i = 0; i < subdomains.size(); i++)
    {
      subdomain_coefs.Append(subdomains[i].scalar_coefficients.Get(scalar_property_name));
    }
    if (!scalars.Has(scalar_property_name))
    {
      scalars.Register(
          scalar_property_name, new mfem::PWCoefficient(subdomain_ids, subdomain_coefs), true);
    }
  }
}
} // namespace hephaestus

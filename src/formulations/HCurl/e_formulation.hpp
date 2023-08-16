#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EFormulation : public hephaestus::HCurlFormulation {
public:
  EFormulation(const std::string &magnetic_reluctivity_name,
               const std::string &magnetic_permeability_name,
               const std::string &electric_conductivity_name,
               const std::string &e_field_name);

  ~EFormulation(){};

  virtual void RegisterCoefficients() override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string &_magnetic_reluctivity_name =
      hephaestus::HCurlFormulation::_alpha_coef_name;
  const std::string &_electric_conductivity_name =
      hephaestus::HCurlFormulation::_beta_coef_name;
};

} // namespace hephaestus

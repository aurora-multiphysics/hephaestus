#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HFormulation : public hephaestus::HCurlFormulation {
public:
  HFormulation(const std::string &electric_resistivity_name,
               const std::string &electric_conductivity_name,
               const std::string &magnetic_permeability_name,
               const std::string &h_field_name);

  ~HFormulation(){};

  virtual void RegisterAuxSolvers() override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _electric_conductivity_name;
  const std::string &_electric_resistivity_name =
      hephaestus::HCurlFormulation::_alpha_coef_name;
  const std::string &_magnetic_permeability_name =
      hephaestus::HCurlFormulation::_beta_coef_name;
};
} // namespace hephaestus

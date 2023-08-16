#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class AFormulation : public hephaestus::HCurlFormulation {
public:
  AFormulation(const std::string &magnetic_reluctivity_name,
               const std::string &magnetic_permeability_name,
               const std::string &electric_conductivity_name,
               const std::string &magnetic_vector_potential_name);

  ~AFormulation(){};

  virtual void RegisterAuxSolvers() override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string &_magnetic_reluctivity_name =
      hephaestus::HCurlFormulation::_alpha_coef_name;
  const std::string &_electric_conductivity_name =
      hephaestus::HCurlFormulation::_beta_coef_name;
};
} // namespace hephaestus

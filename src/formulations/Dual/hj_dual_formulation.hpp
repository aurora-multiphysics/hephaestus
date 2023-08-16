#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HJDualFormulation : public hephaestus::DualFormulation {
public:
  HJDualFormulation(const std::string &electric_resistivity_name,
                    const std::string &electric_conductivity_name,
                    const std::string &magnetic_permeability_name,
                    const std::string &h_field_name,
                    const std::string &j_field_name);

  ~HJDualFormulation(){};

  virtual void RegisterCoefficients() override;

protected:
  const std::string _electric_conductivity_name;
  const std::string &_electric_resistivity_name =
      hephaestus::DualFormulation::_alpha_coef_name;
  const std::string &_magnetic_permeability_name =
      hephaestus::DualFormulation::_beta_coef_name;
};

} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "statics_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class MagnetostaticFormulation : public hephaestus::StaticsFormulation {
public:
  MagnetostaticFormulation(const std::string &magnetic_reluctivity_name,
               const std::string &magnetic_permeability_name,
               const std::string &magnetic_vector_potential_name);

  ~MagnetostaticFormulation(){};

  virtual void RegisterAuxSolvers() override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string &_magnetic_reluctivity_name =
      hephaestus::StaticsFormulation::_alpha_coef_name;
};
} // namespace hephaestus

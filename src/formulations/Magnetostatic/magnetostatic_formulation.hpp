#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"
#include "statics_formulation.hpp"

namespace hephaestus
{

class MagnetostaticFormulation : public hephaestus::StaticsFormulation
{
public:
  MagnetostaticFormulation(const std::string & magnetic_reluctivity_name,
                           const std::string & magnetic_permeability_name,
                           const std::string & magnetic_vector_potential_name);

  ~MagnetostaticFormulation(){};

  // Enable auxiliary calculation of B ∈ H(div)
  virtual void registerMagneticFluxDensityAux(const std::string & b_field_name) override;

  // Enable auxiliary calculation of H ∈ H(curl)
  virtual void registerMagneticFieldAux(const std::string & h_field_name) override;

  // Enable auxiliary calculation of F ∈ L2
  virtual void registerLorentzForceDensityAux(const std::string & f_field_name,
                                              const std::string & b_field_name,
                                              const std::string & j_field_name) override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string & _magnetic_reluctivity_name = hephaestus::StaticsFormulation::_alpha_coef_name;
};
} // namespace hephaestus

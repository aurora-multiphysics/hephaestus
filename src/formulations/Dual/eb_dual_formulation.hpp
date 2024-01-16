#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus
{

class EBDualFormulation : public hephaestus::DualFormulation
{
public:
  EBDualFormulation(const std::string & magnetic_reluctivity_name,
                    const std::string & magnetic_permeability_name,
                    const std::string & electric_conductivity_name,
                    const std::string & e_field_name,
                    const std::string & b_field_name);

  ~EBDualFormulation() override {}

  // Enable auxiliary calculation of J ∈ H(div)
  virtual void registerCurrentDensityAux(const std::string & j_field_name) override;

  // Enable auxiliary calculation of F ∈ L2
  virtual void registerLorentzForceDensityAux(const std::string & f_field_name,
                                              const std::string & b_field_name,
                                              const std::string & j_field_name) override;

  // Enable auxiliary calculation of P ∈ L2
  virtual void registerJouleHeatingDensityAux(const std::string & p_field_name,
                                              const std::string & e_field_name,
                                              const std::string & conductivity_coef_name) override;

  virtual void RegisterCoefficients() override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string & _magnetic_reluctivity_name = hephaestus::DualFormulation::_alpha_coef_name;
  const std::string & _electric_conductivity_name = hephaestus::DualFormulation::_beta_coef_name;
};

} // namespace hephaestus

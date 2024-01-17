#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus
{

class AFormulation : public hephaestus::HCurlFormulation
{
public:
  AFormulation(const std::string & magnetic_reluctivity_name,
               const std::string & magnetic_permeability_name,
               const std::string & electric_conductivity_name,
               const std::string & magnetic_vector_potential_name);

  ~AFormulation() override {}

  // Enable auxiliary calculation of J ∈ H(div)
  virtual void registerCurrentDensityAux(const std::string & j_field_name) override;

  // Enable auxiliary calculation of B ∈ H(div)
  virtual void registerMagneticFluxDensityAux(const std::string & b_field_name) override;

  // Enable auxiliary calculation of E ∈ H(curl)
  virtual void registerElectricFieldAux(const std::string & e_field_name) override;

  // Enable auxiliary calculation of H ∈ H(curl)
  virtual void registerMagneticFieldAux(const std::string & h_field_name) override;

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
  const std::string & _magnetic_reluctivity_name = hephaestus::HCurlFormulation::_alpha_coef_name;
  const std::string & _electric_conductivity_name = hephaestus::HCurlFormulation::_beta_coef_name;
};
} // namespace hephaestus

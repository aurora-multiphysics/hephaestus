#pragma once
#include "mfem.hpp"

namespace hephaestus
{

// Specifies interface specific to EM formulations.
class EMFormulationInterface
{
public:
  EMFormulationInterface() = default;

  // Enable auxiliary calculation of J ∈ H(div)
  virtual void RegisterCurrentDensityAux(const std::string & j_field_name,
                                         const std::string & external_j_field_name = "")
  {
    MFEM_ABORT("Current density auxsolver not available for this formulation");
  }

  // Enable auxiliary calculation of B ∈ H(div)
  virtual void RegisterMagneticFluxDensityAux(const std::string & b_field_name,
                                              const std::string & external_b_field_name = "")
  {
    MFEM_ABORT("Magnetic flux density auxsolver not available for this formulation");
  }

  // Enable auxiliary calculation of E ∈ H(curl)
  virtual void RegisterElectricFieldAux(const std::string & e_field_name,
                                        const std::string & external_e_field_name = "")
  {
    MFEM_ABORT("Electric field auxsolver not available for this formulation");
  }

  // Enable auxiliary calculation of H ∈ H(curl)
  virtual void RegisterMagneticFieldAux(const std::string & h_field_name,
                                        const std::string & external_h_field_name = "")
  {
    MFEM_ABORT("Magnetic field auxsolver not available for this formulation");
  }

  // Enable auxiliary calculation of F ∈ L2
  virtual void RegisterLorentzForceDensityAux(const std::string & f_field_name,
                                              const std::string & b_field_name,
                                              const std::string & j_field_name)
  {
    MFEM_ABORT("Lorentz force auxsolver not available for this formulation");
  }

  // Enable auxiliary calculation of P ∈ L2
  virtual void RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                              const std::string & e_field_name,
                                              const std::string & j_field_name)
  {
    MFEM_ABORT("Joule heating auxsolver not available for this formulation");
  }
};
} // namespace hephaestus

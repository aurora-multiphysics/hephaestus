#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {
/**
 * Solves:
 * ∇×(ν∇×A) + σdA/dt = Jᵉ
 *
 * in weak form
 * (ν∇×A, ∇×A') + (σdA/dt, A') - (Jᵉ, A') - <(ν∇×A)×n, A'>  = 0
 *
 * where:
 * reluctivity ν = 1/μ
 * electrical_conductivity σ=1/ρ
 * Magnetic vector potential A
 * Electric field, E = ρJᵉ -dA/dt
 * Magnetic flux density, B = ∇×A
 * Magnetic field H = ν∇×A
 * Current density J = Jᵉ -σdA/dt
 */
class AFormulation : public hephaestus::HCurlFormulation {
public:
  AFormulation();

  ~AFormulation(){};

  virtual void RegisterAuxSolvers() override;

  virtual void RegisterCoefficients() override;
};
} // namespace hephaestus

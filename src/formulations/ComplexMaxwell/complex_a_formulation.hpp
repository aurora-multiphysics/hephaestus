#pragma once
#include "complex_maxwell_formulation.hpp"

namespace hephaestus {
/*
Formulation for solving:
∇×(ν∇×A) + iωσA - ω²εA = J

via the weak form:
(ν∇×A, ∇×A') + (iωσA, A') - (ω²εA, A') - <(ν∇×A)×n, A'> = (J, A')

where
A ∈ H(curl) is the magnetic vector potential
J ∈ H(div) is the divergence-free source current density
ω is the angular frequency
ν is the reluctivity (reciprocal of magnetic permeability; 1/µ)
σ is the electric conductivity
ε is the electric permittivity

Dirichlet boundaries strongly constrain n×n×A
Integrated boundaries weakly constrain (ν∇×A)×n = H×n
Robin boundaries weakly constrain (ν∇×A)×n + γ(n×n×A) = F

Divergence cleaning (such as via Helmholtz projection)
should be performed on J before use in this operator.
*/
class ComplexAFormulation : public hephaestus::ComplexMaxwellFormulation {

public:
  ComplexAFormulation(const std::string &magnetic_reluctivity_name,
                      const std::string &electric_conductivity_name,
                      const std::string &dielectric_permittivity_name,
                      const std::string &frequency_coef_name,
                      const std::string &magnetic_vector_potential_name);

  virtual void RegisterAuxSolvers() override;
};
} // namespace hephaestus

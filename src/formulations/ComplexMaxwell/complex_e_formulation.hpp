#pragma once
#include "complex_maxwell_formulation.hpp"

namespace hephaestus {
/*
Formulation for solving:
∇×(ν∇×E꜀) + iωσE꜀ - ω²εE꜀ = -iωJ꜀ᵉ

via the weak form:
(ν∇×E꜀, ∇×E꜀') + (iωσE꜀, E꜀') - (ω²εE꜀, E꜀') - <(ν∇×E꜀)×n, E꜀'> = -(iωJ꜀ᵉ, E꜀')

where
E꜀ ∈ H(curl) is the complex electric field
J꜀ᵉ ∈ H(div) is the divergence-free source current density
ω is the angular frequency
ν is the reluctivity (reciprocal of magnetic permeability; 1/µ)
σ is the electric conductivity
ε is the electric permittivity

Dirichlet boundaries strongly constrain n×n×E꜀
Integrated boundaries weakly constrain (ν∇×E꜀)×n = n×dB꜀/dt
Robin boundaries weakly constrain (ν∇×E꜀)×n + γ(n×n×E꜀) = F

Divergence cleaning (such as via Helmholtz projection)
should be performed on J꜀ᵉ before use in this operator.
*/
class ComplexEFormulation : public hephaestus::ComplexMaxwellFormulation {
public:
  ComplexEFormulation();
};
} // namespace hephaestus

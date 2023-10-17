#pragma once
#include "scaled_vector_gridfunction_aux.hpp"

namespace hephaestus {

// Scale the curl of a gridfunction in H(Curl) by a scalar Coefficient, and
// store the result. Suitable for solving for H(Div) or H(Curl) conforming
// fields for expressions like E = ρ∇×H
class ScaledCurlVectorGridFunctionAux : public ScaledVectorGridFunctionAux {
public:
  ScaledCurlVectorGridFunctionAux(const hephaestus::InputParameters &params);
  virtual void buildMixedBilinearForm() override;
};
} // namespace hephaestus

#pragma once
#include "steady_state_problem_builder.hpp"

namespace hephaestus {

// Specifies output interfaces of a frequency-domain EM formulation.
class FrequencyDomainFormulation
    : public hephaestus::SteadyStateProblemBuilder {
protected:
  mfem::ConstantCoefficient *freqCoef;
public:
  FrequencyDomainFormulation();
};
} // namespace hephaestus

#pragma once
#include "frequency_domain_problem_builder.hpp"

namespace hephaestus {

// Specifies output interfaces of a frequency-domain EM formulation.
class FrequencyDomainFormulation
    : public hephaestus::FrequencyDomainProblemBuilder {

public:
  FrequencyDomainFormulation();
};
} // namespace hephaestus

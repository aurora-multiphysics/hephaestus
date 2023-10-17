#pragma once
#include "steady_state_problem_builder.hpp"

namespace hephaestus {

// Specifies output interfaces of a time-independent problem formulation.
class SteadyStateFormulation : public hephaestus::SteadyStateProblemBuilder {
public:
  SteadyStateFormulation();
};
} // namespace hephaestus

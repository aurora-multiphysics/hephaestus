#pragma once
#include "time_domain_problem_builder.hpp"

namespace hephaestus {

// Abstract Factory class of a time-domain EM formulation.
class TransientFormulation : public hephaestus::TransientProblemBuilder {
public:
  TransientFormulation();
};
} // namespace hephaestus

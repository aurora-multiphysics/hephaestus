#pragma once
#include "em_formulation_interface.hpp"
#include "time_domain_problem_builder.hpp"

namespace hephaestus
{

// Abstract Factory class of a time-domain EM formulation.
class TimeDomainEMFormulation : public hephaestus::TimeDomainProblemBuilder,
                                public hephaestus::EMFormulationInterface
{
public:
  TimeDomainEMFormulation();
};
} // namespace hephaestus

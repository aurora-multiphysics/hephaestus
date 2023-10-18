#pragma once
#include "em_formulation_interface.hpp"
#include "steady_state_problem_builder.hpp"

namespace hephaestus {

// Specifies output interfaces of a frequency-domain EM formulation.
class FrequencyDomainEMFormulation
    : public hephaestus::SteadyStateProblemBuilder,
      public hephaestus::EMFormulationInterface {
protected:
  mfem::ConstantCoefficient *freqCoef;

public:
  FrequencyDomainEMFormulation();
};
} // namespace hephaestus

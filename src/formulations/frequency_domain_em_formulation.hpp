#pragma once
#include "complex_em_formulation_interface.hpp"
#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

// Specifies output interfaces of a frequency-domain EM formulation.
class FrequencyDomainEMFormulation : public hephaestus::SteadyStateProblemBuilder,
                                     public hephaestus::ComplexEMFormulationInterface
{
public:
  FrequencyDomainEMFormulation();
  ~FrequencyDomainEMFormulation() override = default;

protected:
  std::shared_ptr<mfem::ConstantCoefficient> _freq_coef{nullptr};
};
} // namespace hephaestus

#include "frequency_domain_problem_builder.hpp"

namespace hephaestus {
FrequencyDomainProblem::FrequencyDomainProblem(
    const hephaestus::InputParameters &params)
    : Problem(params),
      formulation(
          params.GetParam<hephaestus::SteadyFormulation *>("Formulation")) {}

} // namespace hephaestus

#include "time_domain_problem_builder.hpp"

namespace hephaestus {
TransientProblem::TransientProblem(const hephaestus::InputParameters &params)
    : Problem(params),
      formulation(
          params.GetParam<hephaestus::TransientFormulation *>("Formulation")) {}

} // namespace hephaestus

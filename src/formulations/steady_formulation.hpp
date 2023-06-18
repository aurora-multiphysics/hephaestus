#pragma once
#include "frequency_domain_problem_builder.hpp"

namespace hephaestus {

// Specifies output interfaces of a frequency-domain EM formulation.
class SteadyFormulation : public hephaestus::FrequencyDomainProblemBuilder {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;
  std::string frequency_coef_name, permittivity_coef_name,
      reluctivity_coef_name, conductivity_coef_name, h_curl_var_name;

public:
  SteadyFormulation();
};
} // namespace hephaestus

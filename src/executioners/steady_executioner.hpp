#pragma once
#include "executioner_base.hpp"
#include "frequency_domain_problem_builder.hpp"

namespace hephaestus {

class SteadyExecutioner : public Executioner {
public:
  SteadyExecutioner() = default;
  explicit SteadyExecutioner(const hephaestus::InputParameters &params);

  void Init() override;

  void Solve() const override;

  void Execute() const override;

  std::unique_ptr<hephaestus::FrequencyDomainProblem> problem;
};

} // namespace hephaestus

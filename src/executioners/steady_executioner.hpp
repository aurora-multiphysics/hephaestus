#pragma once
#include "executioner_base.hpp"

namespace hephaestus {

class SteadyExecutioner : public Executioner {
private:
  std::unique_ptr<hephaestus::FrequencyDomainOperator> fd_operator;

public:
  SteadyExecutioner() = default;
  explicit SteadyExecutioner(const hephaestus::InputParameters &params);

  void Init() override;

  void Solve() const override;

  void Execute() const override;

  hephaestus::SteadyFormulation *formulation;
};

} // namespace hephaestus

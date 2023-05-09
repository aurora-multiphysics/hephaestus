#pragma once
#include "executioner_base.hpp"

namespace hephaestus {

class SteadyExecutioner : public ExecutionerBase {
private:
public:
  SteadyExecutioner() = default;
  explicit SteadyExecutioner(const hephaestus::InputParameters &params);

  void Init() override;

  void Solve() const override;

  hephaestus::SteadyFormulation *formulation;
};

} // namespace hephaestus

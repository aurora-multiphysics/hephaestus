#pragma once
#include "executioner.hpp"

namespace hephaestus {

class TransientExecutioner : public ExecutionerBase {
private:
  mutable double t_step;  // Time step
  double t_initial;       // Start time
  double t_final;         // End time
  mutable double t;       // Current time
  int vis_steps;          // Number of cyces between each output update
  mutable bool last_step; // Flag to check if current step is final

public:
  TransientExecutioner() = default;
  explicit TransientExecutioner(const hephaestus::InputParameters &params);

  void Init() override;

  void Step(double dt, int it) const;

  void Solve() const override;

  hephaestus::TransientFormulation *formulation;
};

} // namespace hephaestus

#pragma once
#include "executioner_base.hpp"
#include "time_domain_problem_builder.hpp"

namespace hephaestus
{

class TransientExecutioner : public Executioner
{
private:
  double t_initial;       // Start time
  double t_final;         // End time
  mutable double t;       // Current time
  mutable int it;         // Time index
  int vis_steps;          // Number of cyces between each output update
  mutable bool last_step; // Flag to check if current step is final

public:
  mutable double t_step; // Time step

  TransientExecutioner() = default;
  explicit TransientExecutioner(const hephaestus::InputParameters & params);

  void Step(double dt, int it) const;

  void Solve() const override;

  void Execute() const override;

  hephaestus::TimeDomainProblem * problem;
};

} // namespace hephaestus

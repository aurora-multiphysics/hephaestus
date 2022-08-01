#pragma once
#include "auxkernels.hpp"
#include "inputs.hpp"
#include "postprocessors.hpp"
#include "variables.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

class TransientExecutioner {
public:
  TransientExecutioner() {}
  TransientExecutioner(const hephaestus::InputParameters &params);
  void Solve(const hephaestus::InputParameters &params);

  double dt;
  double t_initial;
  double t_final;
  double t;
  int vis_steps;
  bool visualization;
  bool last_step;
};

} // namespace hephaestus

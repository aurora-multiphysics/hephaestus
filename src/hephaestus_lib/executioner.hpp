#pragma once
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

class Executioner {
public:
  Executioner() {}
  Executioner(const std::string &executioner_type, const double time_step,
              const double end_time);

  std::string type;
  double dt;
  double t_final;
};

} // namespace hephaestus
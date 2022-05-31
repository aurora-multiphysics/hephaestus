#pragma once
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

class Executioner {
public:
  Executioner() {}
  Executioner(const std::string &type_, const double dt_,
              const double t_initial_, const double t_final_);

  std::string type;
  double dt;
  double t_initial;
  double t_final;
};

} // namespace hephaestus

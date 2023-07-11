#pragma once
#include "inputs.hpp"

namespace hephaestus {

class Executioner {
protected:
  bool visualization; // Flag to control whether GLVis visualisation is required

public:
  Executioner() = default;
  explicit Executioner(const hephaestus::InputParameters &params);

  // Initialise owned objects
  virtual void Init() = 0;

  // Solve the current system of equations
  virtual void Solve() const = 0;

  // Execute solution strategy including any timestepping
  virtual void Execute() const = 0;

  // Enable output to GLVis
  void EnableVisualisation() { visualization = true; };
};

} // namespace hephaestus

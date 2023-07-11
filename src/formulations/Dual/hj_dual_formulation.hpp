#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HJDualFormulation : public hephaestus::DualFormulation {
public:
  HJDualFormulation();

  ~HJDualFormulation(){};

  virtual void RegisterCoefficients() override;
};

} // namespace hephaestus

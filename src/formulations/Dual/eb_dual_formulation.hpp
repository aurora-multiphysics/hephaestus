#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EBDualFormulation : public hephaestus::DualFormulation {
public:
  EBDualFormulation();

  ~EBDualFormulation(){};

  virtual void RegisterCoefficients() override;
};

} // namespace hephaestus

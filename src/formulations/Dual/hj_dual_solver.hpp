#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HJDualFormulation : public hephaestus::DualFormulation {
public:
  HJDualFormulation();

  ~HJDualFormulation(){};

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};

} // namespace hephaestus

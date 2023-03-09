#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EBDualFormulation : public hephaestus::DualFormulation {
public:
  EBDualFormulation();

  ~EBDualFormulation(){};

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};

} // namespace hephaestus

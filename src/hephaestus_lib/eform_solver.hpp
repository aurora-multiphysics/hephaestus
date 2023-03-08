#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EFormulation : public hephaestus::HCurlFormulation {
public:
  EFormulation();

  ~EFormulation(){};

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};

} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class EFormulation : public hephaestus::HCurlFormulation {
public:
  EFormulation();

  ~EFormulation(){};

  virtual void RegisterCoefficients() override;
};

} // namespace hephaestus

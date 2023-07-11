#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HFormulation : public hephaestus::HCurlFormulation {
public:
  HFormulation();

  ~HFormulation(){};

  virtual void RegisterAuxSolvers() override;

  virtual void RegisterCoefficients() override;
};
} // namespace hephaestus

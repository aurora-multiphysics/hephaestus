#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class AFormulation : public hephaestus::HCurlFormulation {
public:
  AFormulation();

  ~AFormulation(){};

  virtual void RegisterAuxSolvers() override;

  virtual void RegisterCoefficients() override;
};
} // namespace hephaestus

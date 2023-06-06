#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HFormulation : public hephaestus::HCurlFormulation {
public:
  HFormulation();

  ~HFormulation(){};

  virtual void
  RegisterAuxSolvers(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxSolvers &preprocessors) override;

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};
} // namespace hephaestus

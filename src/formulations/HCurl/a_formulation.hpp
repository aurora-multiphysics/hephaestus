#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class AFormulation : public hephaestus::HCurlFormulation {
public:
  AFormulation();

  ~AFormulation(){};

  virtual void
  RegisterAuxSolvers(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxSolvers &auxsolvers) override;

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};
} // namespace hephaestus

#pragma once
#include "div_free_source.hpp"
#include "scalar_potential_source.hpp"

namespace hephaestus {

class Sources : public mfem::NamedFieldsMap<hephaestus::Source> {
public:
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &domain_properties);
  void Apply(mfem::ParLinearForm *lf);
  void SubtractSources(mfem::ParGridFunction *gf);
};

}; // namespace hephaestus

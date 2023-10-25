#pragma once
#include "div_free_source.hpp"
#include "scalar_potential_source.hpp"
#include "closed_coil.hpp"

namespace hephaestus {

class Sources : public mfem::NamedFieldsMap<hephaestus::Source> {
public:
  void Init(hephaestus::GridFunctions &gridfunctions,
            const hephaestus::FESpaces &fespaces, hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &coefficients);
  void Apply(mfem::ParLinearForm *lf);
  void SubtractSources(mfem::ParGridFunction *gf);
};

}; // namespace hephaestus

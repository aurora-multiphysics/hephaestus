#include "sources.hpp"

namespace hephaestus {

void Sources::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map, hephaestus::Coefficients &coefficients) {
  for (const auto &[name, source] : GetMap()) {
    source->Init(variables, fespaces, bc_map, coefficients);
  }
}

void Sources::Apply(mfem::ParLinearForm *lf) {
  for (const auto &[name, source] : GetMap()) {
    source->Apply(lf);
  }
}

void Sources::SubtractSources(mfem::ParGridFunction *gf) {
  for (const auto &[name, source] : GetMap()) {
    source->SubtractSource(gf);
  }
}

} // namespace hephaestus

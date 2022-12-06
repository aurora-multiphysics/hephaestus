#pragma once
#include <any>
#include <fstream>
#include <iostream>
#include <memory>

#include "variables.hpp"

namespace hephaestus {

class Source {
public:
  Source() {}
  virtual void AddSource(mfem::LinearForm lf){};
};

class DivFreeVolumetricSource : public hephaestus::Source {
public:
  DivFreeVolumetricSource(mfem::VectorFunctionCoefficient *src_coef);
  void AddSource(mfem::LinearForm lf);

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field
  mfem::ParBilinearForm *hCurlMass;
  mfem::common::DivergenceFreeProjector *divFreeProj;
};

} // namespace hephaestus

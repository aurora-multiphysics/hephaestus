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
  virtual void
  Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces){};
  virtual void ApplySource(mfem::LinearForm *lf){};
};

class DivFreeVolumetricSource : public hephaestus::Source {
public:
  DivFreeVolumetricSource(const hephaestus::InputParameters &params);
  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces)
      override;
  void ApplySource(mfem::LinearForm *lf) override;

  std::string src_gf_name;
  std::string hcurl_fespace_name;
  std::string h1_fespace_name;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field
  mfem::ParBilinearForm *h_curl_mass;
  mfem::common::DivergenceFreeProjector *divFreeProj;
  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
};

} // namespace hephaestus

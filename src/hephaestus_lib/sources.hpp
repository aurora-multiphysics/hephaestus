#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"
#include "materials.hpp"

namespace hephaestus {

class Source {
public:
  Source() {}
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::DomainProperties &domain_properties){};
  virtual void ApplySource(mfem::LinearForm *lf){};
};

class Sources : public mfem::NamedFieldsMap<hephaestus::Source> {
public:
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::DomainProperties &domain_properties);
  void ApplySources(mfem::LinearForm *lf);
};

class DivFreeVolumetricSource : public hephaestus::Source {
public:
  DivFreeVolumetricSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::DomainProperties &domain_properties) override;
  void ApplySource(mfem::LinearForm *lf) override;

  std::string src_gf_name;
  std::string src_coef_name;
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

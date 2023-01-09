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
  virtual void ApplySource(mfem::ParLinearForm *lf){};
};

class Sources : public mfem::NamedFieldsMap<hephaestus::Source> {
public:
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::DomainProperties &domain_properties);
  void ApplySources(mfem::ParLinearForm *lf);
};

class DivFreeVolumetricSource : public hephaestus::Source {
public:
  DivFreeVolumetricSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::DomainProperties &domain_properties) override;
  void ApplySource(mfem::ParLinearForm *lf) override;

  std::string src_gf_name;
  std::string src_coef_name;
  std::string hcurl_fespace_name;
  std::string h1_fespace_name;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field
  mfem::common::DivergenceFreeProjector *divFreeProj;
  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::Solver *solver;

  const hephaestus::InputParameters solver_options;

protected:
  class DefaultGMRESSolver : public mfem::HypreGMRES {
  public:
    DefaultGMRESSolver(const hephaestus::InputParameters &params,
                       const mfem::HypreParMatrix &M)
        : mfem::HypreGMRES(M), amg(M),
          tol(params.GetOptionalParam<double>("Tolerance", 1e-12)),
          max_iter(params.GetOptionalParam<int>("MaxIter", 200)),
          print_level(params.GetOptionalParam<int>("PrintLevel", 0)) {

      SetPrintLevel(0);
      SetTol(tol);
      SetMaxIter(max_iter);
      SetPrintLevel(print_level);
      SetPreconditioner(amg);
    }
    mfem::HypreBoomerAMG amg;
    double tol;
    int max_iter;
    int print_level;
  };
};
}; // namespace hephaestus

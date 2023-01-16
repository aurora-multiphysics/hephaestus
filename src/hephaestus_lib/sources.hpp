#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "materials.hpp"

namespace hephaestus {

class Source {
public:
  Source() {}
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties){};
  virtual void ApplySource(mfem::ParLinearForm *lf){};
};

class Sources : public mfem::NamedFieldsMap<hephaestus::Source> {
public:
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties);
  void ApplySources(mfem::ParLinearForm *lf);
};

class ScalarPotentialSource : public hephaestus::Source {
public:
  ScalarPotentialSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties) override;
  void ApplySource(mfem::ParLinearForm *lf) override;
  void buildM1(mfem::Coefficient *Sigma);
  void buildGrad();

  std::string src_gf_name;
  std::string potential_gf_name;
  std::string hcurl_fespace_name;
  std::string h1_fespace_name;
  std::string beta_coef_name;
  const hephaestus::InputParameters solver_options;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParGridFunction *p_; // Potential
  hephaestus::BCMap *_bc_map;
  mfem::Coefficient *betaCoef;

  mfem::ParBilinearForm *a0;
  mfem::ParBilinearForm *m1;
  mfem::HypreParMatrix *A0;
  mfem::Vector *X0, *B0;
  mutable mfem::HypreSolver *amg_a0;
  mutable hephaestus::DefaultH1PCGSolver *a0_solver;

  mfem::ParDiscreteLinearOperator *grad;

  mfem::ParLinearForm *b0;
  mfem::ParGridFunction *grad_p_;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field
  mfem::common::DivergenceFreeProjector *divFreeProj;

  mfem::Solver *solver;
};

class DivFreeVolumetricSource : public hephaestus::Source {
public:
  DivFreeVolumetricSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
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
};
}; // namespace hephaestus

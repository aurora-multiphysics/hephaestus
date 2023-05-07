#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "kernels.hpp"
#include "materials.hpp"

namespace hephaestus {

class Source : public hephaestus::Kernel<mfem::ParLinearForm> {
public:
  Source() {}
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties){};
  virtual void Apply(mfem::ParLinearForm *lf) override = 0;
  virtual void SubtractSource(mfem::ParGridFunction *gf) = 0;
};

class Sources : public mfem::NamedFieldsMap<hephaestus::Source> {
public:
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties);
  void Apply(mfem::ParLinearForm *lf);
  void SubtractSources(mfem::ParGridFunction *gf);
};

class DivFreeSource : public hephaestus::Source {
public:
  DivFreeSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;
  void buildH1Diffusion();
  void buildHCurlMass();
  void buildWeakDiv();
  void buildGrad();

  std::string src_gf_name;
  std::string src_coef_name;
  std::string potential_gf_name;
  std::string hcurl_fespace_name;
  std::string h1_fespace_name;
  const hephaestus::InputParameters solver_options;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParGridFunction *q_; // Potential
  hephaestus::BCMap *_bc_map;
  mfem::Coefficient *betaCoef;

  mfem::ParBilinearForm *a0, *s0_;
  mfem::ParBilinearForm *h_curl_mass;
  mfem::ParMixedBilinearForm *weakDiv_;
  mfem::ParDiscreteLinearOperator *grad, *grad_;

  mfem::HypreParMatrix *A0, *S0_;
  mfem::Vector *X0, *B0;
  mutable mfem::HypreSolver *amg_a0;
  mutable hephaestus::DefaultH1PCGSolver *a0_solver;

  mfem::ParLinearForm *gDiv_;
  mfem::ParGridFunction *grad_p_;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field

  mfem::Solver *solver;
  int irOrder, geom;
};

class ScalarPotentialSource : public hephaestus::Source {
public:
  ScalarPotentialSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;
  void buildH1Diffusion(mfem::Coefficient *Sigma);
  void buildM1(mfem::Coefficient *Sigma);
  void buildHCurlMass();
  void buildWeakDiv();
  void buildGrad();

  std::string src_gf_name;
  std::string src_coef_name;
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

  mfem::ParBilinearForm *a0, *s0_;
  mfem::ParBilinearForm *m1;
  mfem::ParBilinearForm *h_curl_mass;
  mfem::ParMixedBilinearForm *weakDiv_;
  mfem::ParDiscreteLinearOperator *grad, *grad_;

  mfem::HypreParMatrix *A0, *S0_;
  mfem::Vector *X0, *B0;
  mutable mfem::HypreSolver *amg_a0;
  mutable hephaestus::DefaultH1PCGSolver *a0_solver;

  mfem::ParLinearForm *b0;
  mfem::ParGridFunction *grad_p_, *xDiv_;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field

  mfem::Solver *solver;
  int irOrder, geom;
};

class DivFreeVolumetricSource : public hephaestus::Source {
public:
  DivFreeVolumetricSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;

  std::string src_gf_name;
  std::string src_coef_name;
  std::string hcurl_fespace_name;
  std::string h1_fespace_name;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *div_free_src_gf; // Source field
  mfem::common::DivergenceFreeProjector *divFreeProj;
  mfem::VectorFEMassIntegrator *h_curl_mass_integ;
  mfem::ParBilinearForm *h_curl_mass;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::Solver *solver;

  const hephaestus::InputParameters solver_options;
};
}; // namespace hephaestus

#pragma once
#include "source.hpp"

namespace hephaestus {

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

}; // namespace hephaestus

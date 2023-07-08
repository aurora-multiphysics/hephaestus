#pragma once
#include "source_base.hpp"

namespace hephaestus {

class DivFreeSource : public hephaestus::Source {
public:
  DivFreeSource(const hephaestus::InputParameters &params);
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &coefficients) override;
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
  bool perform_helmholtz_projection;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParGridFunction *q_, *grad_q; // Potential
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
  mfem::ParGridFunction *g; // H(Curl) projection of user specified source
  mfem::ParGridFunction *div_free_src_gf; // Divergence free projected source
  mfem::Solver *solver;
};

}; // namespace hephaestus

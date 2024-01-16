#pragma once
#include "source_base.hpp"

namespace hephaestus
{

class ScalarPotentialSource : public hephaestus::Source
{
public:
  ScalarPotentialSource(const hephaestus::InputParameters & params);

  // Override virtual Source destructor to prevent leaks.
  ~ScalarPotentialSource() override;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParLinearForm * lf) override;
  void SubtractSource(mfem::ParGridFunction * gf) override;
  void buildH1Diffusion(mfem::Coefficient * Sigma);
  void buildM1(mfem::Coefficient * Sigma);
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

  mfem::ParFiniteElementSpace * H1FESpace_;
  mfem::ParFiniteElementSpace * HCurlFESpace_;
  mfem::ParGridFunction * p_; // Potential
  hephaestus::BCMap * _bc_map;
  mfem::Coefficient * betaCoef;

  std::unique_ptr<mfem::ParBilinearForm> a0{nullptr};
  mfem::ParBilinearForm * s0_;

  std::unique_ptr<mfem::ParBilinearForm> m1{nullptr};

  mfem::ParBilinearForm * h_curl_mass;
  mfem::ParMixedBilinearForm * weakDiv_;

  std::unique_ptr<mfem::ParDiscreteLinearOperator> grad{nullptr};
  mfem::ParDiscreteLinearOperator * grad_;

  std::unique_ptr<mfem::HypreParMatrix> A0{nullptr};
  mfem::HypreParMatrix * S0_;

  std::unique_ptr<mfem::Vector> X0{nullptr};
  std::unique_ptr<mfem::Vector> B0{nullptr};

  mutable mfem::HypreSolver * amg_a0;
  mutable std::unique_ptr<hephaestus::DefaultH1PCGSolver> a0_solver{nullptr};

  std::unique_ptr<mfem::ParLinearForm> b0{nullptr};
  mfem::ParGridFunction *grad_p_, *xDiv_;

  mfem::VectorCoefficient * sourceVecCoef;
  mfem::ParGridFunction * div_free_src_gf; // Source field

  mfem::Solver * solver;
  int irOrder, geom;
};

}; // namespace hephaestus

#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

class HelmholtzProjector {
public:
  HelmholtzProjector(const hephaestus::InputParameters &params);

  void Project(hephaestus::GridFunctions &gridfunctions,
               const hephaestus::FESpaces &fespaces, hephaestus::BCMap &bc_map);

  void setForms();
  void setGrad();
  void setBCs();
  void solveLinearSystem();

private:
  std::string hcurl_fespace_name_;
  std::string h1_fespace_name_;
  std::string gf_grad_name_;
  std::string gf_name_;
  hephaestus::InputParameters solver_options_;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParGridFunction *q_;

  // H(Curl) projection of user specified source
  std::unique_ptr<mfem::ParGridFunction> g;

  mfem::ParGridFunction *div_free_src_gf_; // Divergence free projected source

  std::unique_ptr<mfem::ParLinearForm> gDiv_;
  std::unique_ptr<mfem::ParBilinearForm> a0_;
  std::unique_ptr<mfem::ParMixedBilinearForm> weakDiv_;
  std::unique_ptr<mfem::ParDiscreteLinearOperator> grad_;

  mfem::Array<int> ess_bdr_tdofs_;
  hephaestus::BCMap *bc_map_;
};

} // namespace hephaestus

#pragma once
#include "source_base.hpp"

namespace hephaestus {

class DivFreeSource : public hephaestus::Source {
public:
  DivFreeSource(const hephaestus::InputParameters &params);
  void Init(hephaestus::GridFunctions &gridfunctions,
            const hephaestus::FESpaces &fespaces, hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &coefficients) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;
  void buildHCurlMass();

  std::string src_gf_name;
  std::string src_coef_name;
  std::string potential_gf_name;
  std::string hcurl_fespace_name;
  std::string h1_fespace_name;
  const hephaestus::InputParameters solver_options;
  bool perform_helmholtz_projection;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;

  std::unique_ptr<mfem::ParGridFunction> q_; // Potential

  hephaestus::BCMap *bc_map_;
  hephaestus::GridFunctions *gridfunctions_;
  const hephaestus::FESpaces *fespaces_;

  std::unique_ptr<mfem::ParBilinearForm> h_curl_mass;

  mfem::VectorCoefficient *sourceVecCoef;

  // H(Curl) projection of user specified source
  std::unique_ptr<mfem::ParGridFunction> g;

  // Divergence free projected source
  std::unique_ptr<mfem::ParGridFunction> div_free_src_gf;

  mfem::Solver *solver;
};

}; // namespace hephaestus

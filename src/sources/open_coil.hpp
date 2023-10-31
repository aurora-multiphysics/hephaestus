#pragma once
#include "scalar_potential_source.hpp"
#include "source_base.hpp"

namespace hephaestus {

class OpenCoilSolver : public hephaestus::Source {

public:
  OpenCoilSolver(const hephaestus::InputParameters &params,
                   const std::vector<hephaestus::Subdomain> &coil_dom,
                   const std::pair<int,int> electrodes, const int order);

  ~OpenCoilSolver();

  void Init(hephaestus::GridFunctions &gridfunctions,
            const hephaestus::FESpaces &fespaces, hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &coefficients) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;

  private:

  // Parameters
  int order_;
  std::pair<int, int> elec_attrs_;
  std::vector<hephaestus::Subdomain> coil_domains_;
  mfem::ConstantCoefficient *coef1_;
  mfem::ConstantCoefficient *coef0_;
  mfem::Coefficient *Itotal_;

  // Names
  std::string hcurl_fespace_name_;
  std::string h1_fespace_name_;
  std::string J_gf_name_;
  std::string V_gf_name_;
  std::string I_coef_name_;


  // Parent mesh, FE space, and current
  mfem::ParMesh *mesh_parent_;
  mfem::ParGridFunction *J_parent_;
  mfem::ParGridFunction *V_parent_;
  mfem::ParFiniteElementSpace *HCurlFESpace_parent_;
  mfem::ParFiniteElementSpace *H1FESpace_parent_;



};

} // namespace hephaestus
#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class AVSolver : public mfem::TimeDependentOperator {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
  virtual void SetVariableNames();

public:
  AVSolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
           hephaestus::DomainProperties &domain_properties);

  ~AVSolver(){};

  void Init(mfem::Vector &X);

  void buildA1(mfem::Coefficient *sigma, mfem::Coefficient *dtMuInv);
  void buildM1(mfem::Coefficient *sigma);
  void buildCurl(mfem::Coefficient *muInv);
  void buildGrad();

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt);

  void RegisterOutputFields(mfem::DataCollection *dc_);

  void WriteOutputFields(mfem::DataCollection *dc_, int it = 0);

  virtual void WriteConsoleSummary(double t, int it);

  void InitializeGLVis();

  void DisplayToGLVis();
  mfem::common::H1_ParFESpace *H1FESpace_;
  mfem::common::ND_ParFESpace *HCurlFESpace_;
  mfem::common::RT_ParFESpace *HDivFESpace_;

  mfem::Array<int> true_offsets;

  double ElectricLosses() const;

  std::string u_name, p_name, v_name;
  std::string u_display_name, p_display_name, v_display_name;

protected:
  int myid_;
  int num_procs_;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;
  mfem::ParMesh *pmesh_;

  mfem::ParBilinearForm *a0, *a1, *m1;
  mfem::HypreParMatrix *A0, *A1;
  mfem::Vector *X0, *X1, *B0, *B1;

  mfem::ParDiscreteLinearOperator *grad;
  mfem::ParDiscreteLinearOperator *curl;
  mfem::ParMixedBilinearForm *weakCurl;
  mutable mfem::HypreSolver *amg_a0;
  mutable mfem::HyprePCG *pcg_a0;
  mutable mfem::HypreSolver *ams_a1;
  mutable mfem::HyprePCG *pcg_a1;

  // temporary work vectors
  mfem::ParLinearForm *b0, *b1;
  mfem::ParGridFunction e_, de_; // Electric Field (HCurl)
  mfem::ParGridFunction v_, dv_; // Scalar Potential (H1)
  mfem::ParGridFunction *h_;     // Magnetic Field (HCurl)
  mfem::ParGridFunction b_, db_; // Magnetic Flux (HDiv)

  mfem::ParGridFunction *jr_; // Raw Volumetric Current Density (HCurl)
  mfem::ParGridFunction *bd_; // Dual of B (HCurl)
  mfem::ParGridFunction *jd_; // Dual of J, the rhs vector (HCurl)

  double dt_A1;
  mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;      // Reluctivity Coefficient
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef; // Electric Conductivity Coefficient

  // Sockets used to communicate with GLVis
  std::map<std::string, mfem::socketstream *> socks_;
};
} // namespace hephaestus

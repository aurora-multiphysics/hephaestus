#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class AVSolver : public TransientFormulation {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
  virtual void
  SetSourceCoefficient(hephaestus::DomainProperties &domain_properties);
  virtual void SetVariableNames();

public:
  AVSolver(mfem::ParMesh &pmesh, int order,
           mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
           hephaestus::BCMap &bc_map,
           hephaestus::DomainProperties &domain_properties);

  ~AVSolver(){};

  void Init(mfem::Vector &X) override;

  void buildA1(mfem::Coefficient *sigma, mfem::Coefficient *dtMuInv);
  void buildM1(mfem::Coefficient *sigma);
  void buildCurl(mfem::Coefficient *muInv);
  void buildGrad();
  void buildSource();

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;

  void RegisterOutputFields(mfem::DataCollection *dc_) override;

  void WriteOutputFields(mfem::DataCollection *dc_, int it = 0) override;

  virtual void WriteConsoleSummary(double t, int it) override;

  void InitializeGLVis() override;

  void DisplayToGLVis() override;
  mfem::common::H1_ParFESpace *H1FESpace_;
  mfem::common::ND_ParFESpace *HCurlFESpace_;
  mfem::common::RT_ParFESpace *HDivFESpace_;

  mfem::Array<int> block_trueOffsets;

  double ElectricLosses() const;

  std::string u_name, p_name, v_name, e_name;
  std::string u_display_name, p_display_name, v_display_name, e_display_name;

  mfem::ParGridFunction a_, da_; // Magnetic Vector Potential (HCurl)
  mfem::ParGridFunction v_, dv_; // Scalar Potential (H1)
  mfem::ParGridFunction b_, db_; // Magnetic Flux (HDiv)
  mfem::ParGridFunction e_, de_; // Electric Field (HCurl)

protected:
  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;

  mfem::ParBilinearForm *a0, *a1, *a00, *a11, *m1;
  mfem::ParMixedBilinearForm *a01, *a10;
  mfem::BlockOperator *blockAV;
  mfem::BlockDiagonalPreconditioner *blockAVPr;
  mfem::BlockVector *x, *rhs;
  mfem::BlockVector *trueX, *trueRhs;

  mfem::HypreParMatrix *A0, *A1, *A10, *A01, *blockA;
  mfem::Vector *X0, *X1, *B0, *B1;

  mfem::ParDiscreteLinearOperator *grad;
  mfem::ParDiscreteLinearOperator *curl;
  mfem::ParMixedBilinearForm *weakCurl;
  mfem::ParBilinearForm *curlCurl;

  mutable mfem::HypreSolver *amg_a0;
  mutable mfem::HyprePCG *pcg_a0;
  mutable mfem::HypreSolver *ams_a1;
  mutable mfem::HyprePCG *pcg_a1;

  // temporary work vectors
  mfem::ParLinearForm *b0, *b1;

  mfem::ParGridFunction *jr_; // Raw Volumetric Current Density (HCurl)
  mfem::ParGridFunction *bd_; // Dual of B (HCurl)
  mfem::ParGridFunction *jd_; // Dual of J, the rhs vector (HCurl)

  double dt_A1;
  mfem::ConstantCoefficient dtCoef; // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef, negCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;               // Reluctivity Coefficient
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef,
      *negBetaCoef; // Electric Conductivity Coefficient

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *src_gf, *div_free_src_gf; // Source field
  mfem::ParBilinearForm *hCurlMass;
  mfem::common::DivergenceFreeProjector *divFreeProj;

  // Sockets used to communicate with GLVis
  std::map<std::string, mfem::socketstream *> socks_;
};
} // namespace hephaestus

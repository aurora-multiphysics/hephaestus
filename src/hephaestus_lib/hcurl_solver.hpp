#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HCurlSolver : public TransientFormulation {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
  virtual void
  SetSourceCoefficient(hephaestus::DomainProperties &domain_properties);
  virtual void SetVariableNames();

public:
  HCurlSolver(mfem::ParMesh &pmesh, int order,
              mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
              hephaestus::BCMap &bc_map,
              hephaestus::DomainProperties &domain_properties);

  ~HCurlSolver(){};

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

  double ElectricLosses() const;

  std::string u_name, p_name;
  std::string u_display_name, p_display_name;
  mfem::ParGridFunction u_, du_; // HCurl vector field
  mfem::ParGridFunction p_, dp_; // H1 scalar potential
  mfem::ParGridFunction curl_u_; // HDiv Magnetic Flux Density
  std::map<std::string, mfem::socketstream *> socks_;

protected:
  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;

  mfem::ParBilinearForm *a0, *a1, *m1;
  mfem::HypreParMatrix *A0, *A1;
  mfem::Vector *X0, *X1, *B0, *B1;

  mfem::ParDiscreteLinearOperator *grad;
  mfem::ParDiscreteLinearOperator *curl;
  mfem::ParBilinearForm *curlCurl;
  mutable mfem::HypreSolver *amg_a0;
  mutable mfem::HyprePCG *pcg_a0;
  mutable mfem::HypreAMS *ams_a1;
  mutable mfem::HyprePCG *pcg_a1;

  // temporary work vectors
  mfem::ParLinearForm *b0, *b1;

  double dt_A1;
  mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef;

  mfem::VectorCoefficient *sourceVecCoef;
  mfem::ParGridFunction *src_gf, *div_free_src_gf; // Source field
  mfem::ParBilinearForm *hCurlMass;
  mfem::common::DivergenceFreeProjector *divFreeProj;

  // Sockets used to communicate with GLVis
};
} // namespace hephaestus

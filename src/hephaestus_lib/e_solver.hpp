#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"

namespace hephaestus {

class ESolver : public mfem::TimeDependentOperator {
public:
  ESolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
          hephaestus::DomainProperties &domain_properties);

  ~ESolver(){};

  void Init(mfem::Vector &X);

  void buildM1(mfem::PWCoefficient &sigma);
  void buildCurl(mfem::ConstantCoefficient &muInv);
  void buildGrad();

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt);

  void RegisterOutputFields(mfem::DataCollection *dc_);

  void WriteOutputFields(mfem::DataCollection *dc_, int it = 0);

  void InitializeGLVis();

  void DisplayToGLVis();
  mfem::common::H1_ParFESpace *H1FESpace_;
  mfem::common::ND_ParFESpace *HCurlFESpace_;
  mfem::common::RT_ParFESpace *HDivFESpace_;

  mfem::Array<int> true_offsets;

  double ElectricLosses() const;

private:
  int myid_;
  int num_procs_;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;
  mfem::ParMesh *pmesh_;

  // mfem::ParBilinearForm *curlMuInvCurl_;
  // mfem::ParBilinearForm *hCurlMass_;
  // mfem::ParMixedBilinearForm *sigmaGradH1HCurl_;
  // mfem::ParMixedBilinearForm *hCurlH1Stiff;

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
  mfem::ParGridFunction *b0, *b1;

  mfem::ParGridFunction e_, de_; // Electric Field (HCurl)
  mfem::ParGridFunction v_, dv_; // Scalar Potential (H1)
  mfem::ParGridFunction *h_;     // Magnetic Field (HCurl)
  mfem::ParGridFunction b_, db_; // Magnetic Flux (HDiv)

  mfem::ParGridFunction *jr_; // Raw Volumetric Current Density (HCurl)
  mfem::ParGridFunction *bd_; // Dual of B (HCurl)
  mfem::ParGridFunction *jd_; // Dual of J, the rhs vector (HCurl)

  double dt, mu;
  mfem::ConstantCoefficient muInvCoef, dtMuInvCoef; // Reluctivity Coefficient
  mfem::PWCoefficient sigmaCoef; // Electric Conductivity Coefficient

  // Sockets used to communicate with GLVis
  std::map<std::string, mfem::socketstream *> socks_;
};
} // namespace hephaestus

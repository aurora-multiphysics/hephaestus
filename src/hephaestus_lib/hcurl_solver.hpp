#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class HCurlSolver : public TransientFormulation {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
  virtual void RegisterVariables();

public:
  HCurlSolver(mfem::ParMesh &pmesh, int order,
              mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
              mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
              hephaestus::BCMap &bc_map,
              hephaestus::DomainProperties &domain_properties,
              hephaestus::Kernels &kernels,
              hephaestus::InputParameters &solver_options);

  ~HCurlSolver(){};

  void Init(mfem::Vector &X) override;

  void buildA1(mfem::Coefficient *sigma, mfem::Coefficient *dtMuInv);
  void buildCurl(mfem::Coefficient *muInv);

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

  std::string u_name;
  std::string u_display_name;
  std::string alpha_coef_name, beta_coef_name, dtAlpha_coef_name;

  mfem::ParGridFunction u_, du_; // HCurl vector field
  mfem::ParGridFunction p_;      // H1 scalar potential
  mfem::ParGridFunction curl_u_; // HDiv Magnetic Flux Density
  std::map<std::string, mfem::socketstream *> socks_;

protected:
  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::Kernels _kernels;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;
  hephaestus::InputParameters _solver_options;

  mfem::ParBilinearForm *a1;
  mfem::HypreParMatrix *A1;
  mfem::Vector *X1, *B1;

  mfem::ParDiscreteLinearOperator *curl;
  mfem::ParBilinearForm *curlCurl;
  mutable hephaestus::DefaultHCurlPCGSolver *a1_solver;

  // temporary work vectors
  mfem::ParLinearForm *b1;

  double dt_A1;
  mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef;

  // Sockets used to communicate with GLVis
};
} // namespace hephaestus

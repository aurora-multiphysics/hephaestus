#pragma once
#include "../common/pfem_extras.hpp"
#include "steady_formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class LinearElasticSolver : public SteadyFormulation {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
  virtual void RegisterVariables();

public:
  LinearElasticSolver(mfem::ParMesh &pmesh, int order,
              mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
              mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
              hephaestus::BCMap &bc_map,
              hephaestus::DomainProperties &domain_properties,
              hephaestus::Sources &sources,
              hephaestus::InputParameters &solver_options);

  ~LinearElasticSolver(){};

  void Init(mfem::Vector &X) override;

  void buildA1();
  void buildM1(mfem::Coefficient *sigma);
  void buildSource();

  void RegisterOutputFields(mfem::DataCollection *dc_) override;

  void WriteOutputFields(mfem::DataCollection *dc_, int it = 0) override;

  virtual void WriteConsoleSummary(double t, int it) override;

  virtual void NotMult(const mfem::Vector& x, mfem::Vector& y);

  virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const override {};

  void InitializeGLVis() override;

  void DisplayToGLVis() override;
  mfem::common::H1_ParFESpace *H1FESpace_;
  mfem::common::ND_ParFESpace *HCurlFESpace_;
  mfem::common::RT_ParFESpace *HDivFESpace_;

  double ElectricLosses() const;

  std::string u_name, p_name;
  std::string u_display_name, p_display_name;
  mfem::ParGridFunction u_, u_sol_; // HCurl vector field
  mfem::ParGridFunction p_, dp_; // H1 scalar potential
  mfem::ParGridFunction curl_u_; // HDiv Magnetic Flux Density
  std::map<std::string, mfem::socketstream *> socks_;

protected:
  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::Sources &_sources;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;
  hephaestus::InputParameters _solver_options;

  mfem::ParBilinearForm *a0, *a1, *m1;
  mfem::HypreParMatrix *A0, *A1;
  mfem::Vector *X0, *X1, *B0, *B1;

  mutable mfem::HypreSolver *amg_a0;
  mutable mfem::HyprePCG *pcg_a0;
  mutable hephaestus::DefaultH1PCGSolver *a0_solver;
  mutable hephaestus::DefaultHCurlPCGSolver *a1_solver;

  // temporary work vectors
  mfem::ParLinearForm *b0, *b1;

  double dt_A1;
  mfem::PWCoefficient lambda_func; // Lame coefficient
  mfem::PWCoefficient shear_modulus_func; // Lame oefficient
  mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient

  double t = 1;

  // Sockets used to communicate with GLVis
};
} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class AVSolver : public TransientFormulation {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);

public:
  AVSolver(mfem::ParMesh &pmesh, int order,
           mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
           mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
           hephaestus::BCMap &bc_map,
           hephaestus::DomainProperties &domain_properties,
           hephaestus::Sources &sources,
           hephaestus::InputParameters &solver_options);

  ~AVSolver(){};

  void Init(mfem::Vector &X) override;
  virtual void RegisterVariables();
  virtual void RegisterMissingVariables();
  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;

  std::string alpha_coef_name, beta_coef_name;

  mfem::ParGridFunction u_, du_; // HCurl vector field
  mfem::ParGridFunction p_, dp_; // H1 scalar potential
  mfem::ParGridFunction e_;      // HCurl Electric Field
  mfem::ParGridFunction b_;      // HDiv Magnetic Flux Density

protected:
  int myid_;
  int num_procs_;
  const int _order;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::Sources &_sources;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;
  hephaestus::InputParameters _solver_options;

  hephaestus::AVEquationSystem *_equation_system;
  mutable hephaestus::DefaultGMRESSolver *solver;

  mfem::ParBilinearForm *a0, *a1, *m1;
  mfem::ParMixedBilinearForm *a01, *a10;
  mfem::HypreParMatrix *A0, *A1, *A10, *A01, *blockA;
  mfem::Vector *X0, *X1, *B0, *B1;
  mfem::Array<int> block_trueOffsets;

  mfem::ParDiscreteLinearOperator *grad;
  mfem::ParDiscreteLinearOperator *curl;
  mfem::ParBilinearForm *curlCurl;
  mutable mfem::HypreSolver *amg_a0;
  mutable mfem::HyprePCG *pcg_a0;
  mutable mfem::HypreAMS *ams_a0;
  mutable mfem::HyprePCG *pcg_a1;

  // temporary work vectors
  mfem::ParLinearForm *b0, *b1, *b01, *b10;

  double dt_A0, dt_A1;
  mfem::ConstantCoefficient dtCoef; // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef, negCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef, *negBetaCoef;

  // Sockets used to communicate with GLVis
};
} // namespace hephaestus

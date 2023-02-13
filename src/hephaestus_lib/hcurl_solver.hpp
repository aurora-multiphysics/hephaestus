#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class HCurlSolver : public TransientFormulation {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);

public:
  HCurlSolver(mfem::ParMesh &pmesh, int order,
              mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
              mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
              hephaestus::BCMap &bc_map,
              hephaestus::DomainProperties &domain_properties,
              hephaestus::Sources &sources,
              hephaestus::InputParameters &solver_options);

  ~HCurlSolver(){};

  void Init(mfem::Vector &X) override;
  virtual void RegisterVariables();
  virtual void RegisterMissingVariables();
  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;

  std::string u_name;
  mfem::ParGridFunction *u_, *du_; // HCurl vector field
  std::map<std::string, mfem::socketstream *> socks_;

protected:
  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  int _order;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::Sources &_sources;
  hephaestus::BCMap _bc_map;
  hephaestus::DomainProperties _domain_properties;
  hephaestus::InputParameters _solver_options;

  hephaestus::CurlCurlWeakForm *_weak_form;
  mfem::HypreParMatrix *A1;
  mfem::Vector *X1, *B1;

  mutable hephaestus::DefaultHCurlPCGSolver *a1_solver;

  // temporary work vectors
  mfem::ParLinearForm *b1;

  double dt_A1;
  mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef;
};
} // namespace hephaestus

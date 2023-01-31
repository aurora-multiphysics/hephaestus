#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class WeakForm {
public:
  WeakForm(const std::string test_var_name,
           mfem::ParGridFunction &test_variable);

  ~WeakForm(){};

  // All should share the test variable FESpace
  std::string _test_var_name; // TODO: should remove by updating bcmap
  mfem::ParFiniteElementSpace *test_pfes;
  mfem::ParBilinearForm *blf;
  mfem::ParLinearForm *lf;
  mfem::ParNonlinearForm *nlf;
  mfem::NamedFieldsMap<mfem::ParMixedBilinearForm>
      mblfs; // named according to trial variable

  virtual void buildLinearForm(hephaestus::BCMap &bc_map,
                               hephaestus::Sources &sources) = 0;
  virtual void buildBilinearForm() = 0;
  virtual void buildWeakForm(hephaestus::BCMap &bc_map,
                             hephaestus::Sources &sources) = 0;

  virtual void applyBoundaryConditions(hephaestus::BCMap &bc_map);

  virtual void FormLinearSystem(mfem::HypreParMatrix &A, mfem::Vector &X,
                                mfem::Vector &B);
  virtual void RecoverFEMSolution(mfem::Vector &X,
                                  mfem::ParGridFunction &test_variable);

protected:
  // Variables for setting Dirichlet BCs
  mfem::Array<int> ess_tdof_list;
  mfem::ParGridFunction x;
};

/*
Class to store weak form components for time dependent PDEs
*/
class TimeDependentWeakForm : public WeakForm {
public:
  TimeDependentWeakForm(const std::string test_var_name,
                        mfem::ParGridFunction &test_variable);
  ~TimeDependentWeakForm(){};
  virtual void setTimeStep(double dt);
  virtual void updateWeakForm(hephaestus::BCMap &bc_map,
                              hephaestus::Sources &sources);
  mfem::ConstantCoefficient dtCoef; // Coefficient for timestep scaling
};

/*
Class to store weak form components for the equation
∇×(α∇×u) + βdu/dt = s0
where
s0 ∈ H(div) source field
u ∈ H(curl)
with corresponding weak form
(α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
- (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
blf(du/dt_{n+1}, u') = lf(u')
blf(u, u') = (βu, u') + (αdt∇×u, ∇×u')
lf(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>
using
u_{n+1} = u_{n} + dt du/dt_{n+1}
*/
class CurlCurlWeakForm : public TimeDependentWeakForm {
public:
  CurlCurlWeakForm(const std::string test_var_name,
                   mfem::ParGridFunction &test_variable,
                   mfem::ParGridFunction &coupled_variable,
                   mfem::Coefficient *alphaCoef_, mfem::Coefficient *betaCoef_);
  ~CurlCurlWeakForm(){};
  virtual void buildLinearForm(hephaestus::BCMap &bc_map,
                               hephaestus::Sources &sources) override;
  virtual void buildBilinearForm() override;
  virtual void buildWeakForm(hephaestus::BCMap &bc_map,
                             hephaestus::Sources &sources) override;
  virtual void setTimeStep(double dt) override;
  virtual void updateWeakForm(hephaestus::BCMap &bc_map,
                              hephaestus::Sources &sources) override;

  mfem::ParGridFunction &u_;
  mfem::Coefficient *alphaCoef;
  mfem::Coefficient *betaCoef;
  mfem::ParBilinearForm *curlCurl;
  mfem::Coefficient *dtAlphaCoef;
};

// class EquationSystem {
//   virtual void
//   SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
//   virtual void RegisterVariables();

// public:
//   EquationSystem(mfem::ParMesh &pmesh, int order,
//                  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//                  mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//                  hephaestus::BCMap &bc_map,
//                  hephaestus::DomainProperties &domain_properties,
//                  hephaestus::Sources &sources,
//                  hephaestus::InputParameters &solver_options);

//   ~EquationSystem(){};

//   mfem::NamedFieldsMap<hephaestus::WeakForm> eqns; // Named by test variable

//   mfem::HypreParMatrix *A1;
//   mfem::Vector *X1, *B1;

//   void Init(mfem::Vector &X) override;
// };

// class EquationSystem {
//   virtual void
//   SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
//   virtual void RegisterVariables();

// public:
//   EquationSystem(mfem::ParMesh &pmesh, int order,
//                  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//                  mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//                  hephaestus::BCMap &bc_map,
//                  hephaestus::DomainProperties &domain_properties,
//                  hephaestus::Sources &sources,
//                  hephaestus::InputParameters &solver_options);

//   ~EquationSystem(){};

//   mfem::NamedFieldsMap<hephaestus::WeakForm> eqns; // Named by test variable

//   mfem::HypreParMatrix *A1;
//   mfem::Vector *X1, *B1;

//   void Init(mfem::Vector &X) override;

// protected:
//   int myid_;
//   int num_procs_;
//   mfem::ParMesh *pmesh_;
//   mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
//   mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
//   hephaestus::Sources &_sources;
//   hephaestus::BCMap _bc_map;
//   hephaestus::DomainProperties _domain_properties;
//   hephaestus::InputParameters _solver_options;

//   mfem::ParBilinearForm *a1;
//   mfem::HypreParMatrix *A1;
//   mfem::Vector *X1, *B1;

//   mfem::ParDiscreteLinearOperator *curl;
//   mfem::ParBilinearForm *curlCurl;
//   mutable hephaestus::DefaultHCurlPCGSolver *a1_solver;

//   // temporary work vectors
//   mfem::ParLinearForm *b1;

//   double dt_A1;
//   mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
//   mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient
//   mfem::Coefficient *alphaCoef;
//   mfem::Coefficient *dtAlphaCoef;
//   mfem::Coefficient *betaCoef;

//   // Sockets used to communicate with GLVis
// };
} // namespace hephaestus

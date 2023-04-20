#pragma once
#include "../common/pfem_extras.hpp"
#include "auxkernels.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {
// Physical Constants
// Permittivity of Free Space (units F/m)
static const double epsilon0_ = 8.8541878176e-12;
// Permeability of Free Space (units H/m)
static const double mu0_ = 4.0e-7 * M_PI;
static const double freq_ = 9.3e9; // 10/2pi

// Specifies output interfaces of a time-domain EM formulation.
class HertzOperator : public mfem::Operator {
public:
  HertzOperator(mfem::ParMesh &pmesh, int order,
                mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
                mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                hephaestus::BCMap &bc_map,
                hephaestus::DomainProperties &domain_properties,
                hephaestus::Sources &sources,
                hephaestus::InputParameters &solver_options)
      : myid_(0), num_procs_(1), _order(order), pmesh_(&pmesh),
        _fespaces(fespaces), _variables(variables), _bc_map(bc_map),
        _sources(sources), _domain_properties(domain_properties),
        _solver_options(solver_options){};

  ~HertzOperator(){};

  double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
  double port_width_vector[3] = {0.0, 0.0, 10.16e-3};

  mfem::Vector a3Vec;
  mfem::Vector a2xa3;
  mfem::Vector a3xa1;
  double kc;
  double k0;
  std::complex<double> k_;

  virtual void SetVariables();
  virtual void Init(mfem::Vector &X);
  virtual void Solve(mfem::Vector &X);
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override{};

  mfem::Array<int> true_offsets, block_trueOffsets;
  // Vector of names of state variables used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> state_var_names;
  // Vector of names of recognised auxiliary variables that can be calculated
  // from formulation,
  std::vector<std::string> aux_var_names;
  // Vector of names of active auxiliary variables that are being calculated
  // in formulation,
  std::vector<std::string> active_aux_var_names;

  hephaestus::RobinBC *robin_E_bc;
  std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;
  mfem::ParComplexGridFunction *e_;
  mfem::ParComplexLinearForm *jd_;
  mfem::ParSesquilinearForm *a1_;
  mfem::ComplexOperator::Convention conv_ = mfem::ComplexOperator::HERMITIAN;

  int myid_;
  int num_procs_;
  const int _order;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::BCMap &_bc_map;
  hephaestus::Sources &_sources;
  hephaestus::DomainProperties &_domain_properties;
  hephaestus::InputParameters &_solver_options;

  mutable hephaestus::DefaultGMRESSolver *solver = NULL;
  mutable hephaestus::DefaultHCurlPCGSolver *a1_solver = NULL;

  mfem::OperatorHandle blockA;
  mfem::BlockVector trueX, trueRhs;

  mfem::Coefficient *epsCoef_;    // Dielectric Material Coefficient
  mfem::Coefficient *muInvCoef_;  // Dia/Paramagnetic Material Coefficient
  mfem::Coefficient *sigmaCoef_;  // Electrical Conductivity Coefficient
  mfem::Coefficient *etaInvCoef_; // Admittance Coefficient

  mfem::Coefficient *omegaCoef_;     // omega expressed as a Coefficient
  mfem::Coefficient *negOmegaCoef_;  // -omega expressed as a Coefficient
  mfem::Coefficient *omega2Coef_;    // omega^2 expressed as a Coefficient
  mfem::Coefficient *negOmega2Coef_; // -omega^2 expressed as a Coefficient
  mfem::Coefficient *massCoef_;      // -omega^2 epsilon
  mfem::Coefficient *posMassCoef_;   // omega^2 epsilon
  mfem::Coefficient *lossCoef_;      // -omega sigma
  mfem::Coefficient *abcCoef_;       // -omega eta^{-1}
  mfem::Coefficient *posAbcCoef_;    // omega eta^{-1}

  mfem::VectorCoefficient *jrCoef_; // Volume Current Density Function
  mfem::VectorCoefficient *jiCoef_; // Volume Current Density Function
  mfem::VectorCoefficient *erCoef_; // Electric Field Boundary Condition
  mfem::VectorCoefficient *eiCoef_; // Electric Field Boundary Condition

  mfem::Array<int> abcs;
  mfem::Array<int> dbcs;

  // Array of 0's and 1's marking the location of absorbing surfaces
  mfem::Array<int> abc_marker_;

  // Array of 0's and 1's marking the location of Dirichlet boundaries
  mfem::Array<int> dbc_marker_;

  mfem::Array<int> ess_bdr_;
  mfem::Array<int> ess_bdr_tdofs_;
  mfem::Array<int> non_k_bdr_;
};

//
// Specifies output interfaces of a time-domain EM formulation.
class HertzFormulation {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;

public:
  hephaestus::HertzOperator *fd_operator;
  mfem::ConstantCoefficient oneCoef;

  HertzFormulation();

  // virtual hephaestus::TimeDependentEquationSystem *CreateEquationSystem();

  virtual hephaestus::HertzOperator *CreateFrequencyDomainOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options);

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables);

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels){};

  virtual void
  RegisterCoefficients(hephaestus::DomainProperties &domain_properties);
};
} // namespace hephaestus

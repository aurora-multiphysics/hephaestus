// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#pragma once

#ifndef MFEM_HERTZ_SOLVER
#define MFEM_HERTZ_SOLVER
#include "../common/pfem_extras.hpp"
#include "boundary_conditions.hpp"

#ifdef MFEM_USE_MPI

#include <map>
#include <string>

using hephaestus::BCMap;

namespace mfem {

using common::DivergenceFreeProjector;
using common::H1_ParFESpace;
using common::ND_ParFESpace;
using common::ParDiscreteCurlOperator;
using common::ParDiscreteGradOperator;
using common::RT_ParFESpace;

namespace electromagnetics {

// Physical Constants
// Permittivity of Free Space (units F/m)
static const double epsilon0_ = 8.8541878176e-12;
// Permeability of Free Space (units H/m)
static const double mu0_ = 4.0e-7 * M_PI;

// Solver options
struct SolverOptions {
  int maxIter;
  int kDim;
  int printLvl;
  double relTol;

  // Euclid Options
  int euLvl;
};

// class SurfaceCurrent;
class HertzSolver {
public:
  enum PrecondType {
    INVALID_PC = -1,
    DIAG_SCALE = 1,
    PARASAILS = 2,
    EUCLID = 3,
    AMS = 4
  };

  enum SolverType {
    INVALID = -1,
    GMRES = 1,
    FGMRES = 2,
    MINRES = 3,
    SUPERLU = 4,
    STRUMPACK = 5,
    MUMPS = 6
  };

  HertzSolver(ParMesh &pmesh, int order, double freq, HertzSolver::SolverType s,
              SolverOptions &sOpts, HertzSolver::PrecondType p,
              ComplexOperator::Convention conv, Coefficient &epsCoef,
              Coefficient &muInvCoef, Coefficient *sigmaCoef,
              Coefficient *etaInvCoef, hephaestus::BCMap bc_map,
              Array<int> &abcs, Array<int> &dbcs,
              void (*j_r_src)(const Vector &, Vector &),
              void (*j_i_src)(const Vector &, Vector &));
  ~HertzSolver();

  HYPRE_Int GetProblemSize();

  void PrintSizes();

  void Assemble();

  void Update();

  void Solve();

  void GetErrorEstimates(Vector &errors);

  void RegisterVisItFields(VisItDataCollection &visit_dc);

  void WriteVisItFields(int it = 0);

  void InitializeGLVis();

  void DisplayToGLVis();

  void DisplayAnimationToGLVis();

private:
  int myid_;
  int num_procs_;
  int order_;
  int logging_;

  SolverType sol_;
  SolverOptions &solOpts_;
  PrecondType prec_;

  ComplexOperator::Convention conv_;

  bool ownsEtaInv_;

  double freq_;

  ParMesh *pmesh_;

  ND_ParFESpace *HCurlFESpace_;

  Array<int> blockTrueOffsets_;

  ParSesquilinearForm *a1_;
  ParBilinearForm *b1_;

  ParComplexGridFunction *e_; // Complex electric field (HCurl)
  ParGridFunction *e_t_;      // Real electric field (HCurl)
  ParComplexGridFunction *j_; // Complex current density (HCurl)

  ParComplexLinearForm *jd_; // Dual of complex current density (HCurl)

  hephaestus::BCMap bc_map_; //

  Coefficient *epsCoef_;    // Dielectric Material Coefficient
  Coefficient *muInvCoef_;  // Dia/Paramagnetic Material Coefficient
  Coefficient *sigmaCoef_;  // Electrical Conductivity Coefficient
  Coefficient *etaInvCoef_; // Admittance Coefficient

  Coefficient *omegaCoef_;     // omega expressed as a Coefficient
  Coefficient *negOmegaCoef_;  // -omega expressed as a Coefficient
  Coefficient *omega2Coef_;    // omega^2 expressed as a Coefficient
  Coefficient *negOmega2Coef_; // -omega^2 expressed as a Coefficient
  Coefficient *massCoef_;      // -omega^2 epsilon
  Coefficient *posMassCoef_;   // omega^2 epsilon
  Coefficient *lossCoef_;      // -omega sigma
  Coefficient *abcCoef_;       // -omega eta^{-1}
  Coefficient *posAbcCoef_;    // omega eta^{-1}

  VectorCoefficient *jrCoef_; // Volume Current Density Function
  VectorCoefficient *jiCoef_; // Volume Current Density Function
  VectorCoefficient *erCoef_; // Electric Field Boundary Condition
  VectorCoefficient *eiCoef_; // Electric Field Boundary Condition

  void (*j_r_src_)(const Vector &, Vector &);
  void (*j_i_src_)(const Vector &, Vector &);

  // Array of 0's and 1's marking the location of absorbing surfaces
  Array<int> abc_marker_;

  // Array of 0's and 1's marking the location of Dirichlet boundaries
  Array<int> dbc_marker_;

  Array<int> *dbcs_;
  Array<int> ess_bdr_;
  Array<int> ess_bdr_tdofs_;
  Array<int> non_k_bdr_;

  VisItDataCollection *visit_dc_;

  std::map<std::string, socketstream *> socks_;
};

} // namespace electromagnetics

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_HERTZ_SOLVER

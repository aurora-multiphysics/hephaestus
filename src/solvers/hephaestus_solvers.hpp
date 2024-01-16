#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"

namespace hephaestus
{

class DefaultH1PCGSolver : public mfem::HyprePCG
{
public:
  DefaultH1PCGSolver(const hephaestus::InputParameters & params, const mfem::HypreParMatrix & M)
    : mfem::HyprePCG(M),
      amg(M),
      tol(params.GetOptionalParam<float>("Tolerance", 1.0e-9)),
      abstol(params.GetOptionalParam<float>("AbsTolerance", 1e-16)),
      max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
      print_level(params.GetOptionalParam<int>("PrintLevel", 2))
  {

    amg.SetPrintLevel(print_level);
    SetTol(tol);
    SetAbsTol(abstol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
    SetPreconditioner(amg);
  }
  mfem::HypreBoomerAMG amg;
  double tol;
  double abstol;
  int max_iter;
  int print_level;
};

class DefaultJacobiPCGSolver : public mfem::HyprePCG
{
public:
  DefaultJacobiPCGSolver(const hephaestus::InputParameters & params, const mfem::HypreParMatrix & M)
    : mfem::HyprePCG(M),
      jacobi(M),
      tol(params.GetOptionalParam<float>("Tolerance", 1.0e-9)),
      abstol(params.GetOptionalParam<float>("AbsTolerance", 1e-16)),
      max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
      print_level(params.GetOptionalParam<int>("PrintLevel", 2))
  {

    SetTol(tol);
    SetAbsTol(abstol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
    SetPreconditioner(jacobi);
  }
  mfem::HypreDiagScale jacobi;
  double tol;
  double abstol;
  int max_iter;
  int print_level;
};

class DefaultHCurlPCGSolver : public mfem::HyprePCG
{
public:
  DefaultHCurlPCGSolver(const hephaestus::InputParameters & params,
                        const mfem::HypreParMatrix & M,
                        mfem::ParFiniteElementSpace * edge_fespace)
    : mfem::HyprePCG(M),
      ams(M, edge_fespace),
      tol(params.GetOptionalParam<float>("Tolerance", 1.0e-16)),
      abstol(params.GetOptionalParam<float>("AbsTolerance", 1e-16)),
      max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
      print_level(params.GetOptionalParam<int>("PrintLevel", -1))
  {

    ams.SetSingularProblem();
    ams.SetPrintLevel(print_level);
    SetTol(tol);
    SetAbsTol(abstol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
    SetPreconditioner(ams);
  }
  mfem::HypreAMS ams;
  double tol;
  double abstol;
  int max_iter;
  int print_level;
};

class DefaultHCurlFGMRESSolver : public mfem::HypreFGMRES
{
public:
  DefaultHCurlFGMRESSolver(const hephaestus::InputParameters & params,
                           const mfem::HypreParMatrix & M,
                           mfem::ParFiniteElementSpace * edge_fespace)
    : mfem::HypreFGMRES(M),
      ams(M, edge_fespace),
      tol(params.GetOptionalParam<float>("Tolerance", 1e-16)),
      max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 100)),
      k_dim(params.GetOptionalParam<unsigned int>("KDim", 10)),
      print_level(params.GetOptionalParam<int>("PrintLevel", -1))
  {

    ams.SetSingularProblem();
    ams.SetPrintLevel(print_level);
    SetTol(tol);
    SetMaxIter(max_iter);
    SetKDim(k_dim);
    SetPrintLevel(print_level);
    SetPreconditioner(ams);
  }

  DefaultHCurlFGMRESSolver(const hephaestus::InputParameters & params,
                           const mfem::HypreParMatrix & M,
                           mfem::ParFiniteElementSpace * edge_fespace,
                           mfem::HypreParVector interior_nodes)
    : mfem::HypreFGMRES(M),
      ams(M, edge_fespace),
      tol(params.GetOptionalParam<float>("Tolerance", 1e-16)),
      max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 100)),
      k_dim(params.GetOptionalParam<unsigned int>("KDim", 10)),
      print_level(params.GetOptionalParam<int>("PrintLevel", -1))
  {
    HYPRE_Solver solver = ams;
    HYPRE_AMSSetInteriorNodes(solver, interior_nodes);
    ams.SetPrintLevel(print_level);
    SetTol(tol);
    SetMaxIter(max_iter);
    SetKDim(k_dim);
    SetPrintLevel(print_level);
    SetPreconditioner(ams);
  }

  mfem::HypreAMS ams;
  double tol;
  int max_iter;
  int k_dim;
  int print_level;
};

class DefaultGMRESSolver : public mfem::HypreGMRES
{
public:
  DefaultGMRESSolver(const hephaestus::InputParameters & params, const mfem::HypreParMatrix & M)
    : mfem::HypreGMRES(M),
      amg(M),
      tol(params.GetOptionalParam<float>("Tolerance", 1e-16)),
      abstol(params.GetOptionalParam<float>("AbsTolerance", 1e-16)),
      max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
      print_level(params.GetOptionalParam<int>("PrintLevel", -1))
  {

    amg.SetPrintLevel(print_level);
    SetTol(tol);
    SetAbsTol(abstol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
    SetPreconditioner(amg);
  }
  mfem::HypreBoomerAMG amg;
  double tol;
  double abstol;
  int max_iter;
  int print_level;
};

} // namespace hephaestus

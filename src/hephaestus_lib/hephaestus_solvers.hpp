#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"

namespace hephaestus {

class DefaultPCGSolver : public mfem::HyprePCG {
public:
  DefaultPCGSolver(const hephaestus::InputParameters &params,
                   const mfem::HypreParMatrix &M,
                   mfem::ParFiniteElementSpace *edge_fespace)
      : mfem::HyprePCG(M), ams(M, edge_fespace),
        tol(params.GetOptionalParam<double>("Tolerance", 1.0e-16)),
        max_iter(params.GetOptionalParam<int>("MaxIter", 1000)),
        print_level(params.GetOptionalParam<int>("PrintLevel", 0)) {

    ams.SetSingularProblem();
    SetPrintLevel(0);
    SetTol(tol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
    SetPreconditioner(ams);
  }
  mfem::HypreAMS ams;
  double tol;
  int max_iter;
  int print_level;
};

class DefaultGMRESSolver : public mfem::HypreGMRES {
public:
  DefaultGMRESSolver(const hephaestus::InputParameters &params,
                     const mfem::HypreParMatrix &M)
      : mfem::HypreGMRES(M), amg(M),
        tol(params.GetOptionalParam<double>("Tolerance", 1e-12)),
        max_iter(params.GetOptionalParam<int>("MaxIter", 200)),
        print_level(params.GetOptionalParam<int>("PrintLevel", 0)) {

    SetPrintLevel(0);
    SetTol(tol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
    SetPreconditioner(amg);
  }
  mfem::HypreBoomerAMG amg;
  double tol;
  int max_iter;
  int print_level;
};

} // namespace hephaestus

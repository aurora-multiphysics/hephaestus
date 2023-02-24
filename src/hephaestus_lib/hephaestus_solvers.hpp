#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"

namespace hephaestus {

class DefaultH1PCGSolver : public mfem::HyprePCG {
public:
  DefaultH1PCGSolver(const hephaestus::InputParameters &params,
                     const mfem::HypreParMatrix &M)
      : mfem::HyprePCG(M), amg(M),
        tol(params.GetOptionalParam<float>("Tolerance", 1.0e-9)),
        max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
        print_level(params.GetOptionalParam<int>("PrintLevel", -1)) {

    amg.SetPrintLevel(print_level);
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

class DefaultHCurlPCGSolver : public mfem::HyprePCG {
public:
  DefaultHCurlPCGSolver(const hephaestus::InputParameters &params,
                        const mfem::HypreParMatrix &M,
                        mfem::ParFiniteElementSpace *edge_fespace)
      : mfem::HyprePCG(M), ams(M, edge_fespace),
        tol(params.GetOptionalParam<float>("Tolerance", 1.0e-16)),
        max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
        print_level(params.GetOptionalParam<int>("PrintLevel", -1)) {

    ams.SetSingularProblem();
    ams.SetPrintLevel(print_level);
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

class DefaultGMRESSolver : public mfem::GMRESSolver {
public:
  DefaultGMRESSolver(const hephaestus::InputParameters &params,
                     const mfem::HypreParMatrix &M)
      : mfem::GMRESSolver(M.GetComm()),
        atol(params.GetOptionalParam<float>("AbsoluteTolerance", 1.0e-16)),
        max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
        print_level(params.GetOptionalParam<int>("PrintLevel", -1)) {

    SetOperator(M);
    SetAbsTol(atol);
    SetMaxIter(max_iter);
    SetPrintLevel(print_level);
  }
  double atol;
  int max_iter;
  int print_level;
};

// class DefaultGMRESSolver : public mfem::HypreGMRES {
// public:
//   DefaultGMRESSolver(const hephaestus::InputParameters &params,
//                      const mfem::HypreParMatrix &M)
//       : mfem::HypreGMRES(M), amg(M),
//         tol(params.GetOptionalParam<float>("Tolerance", 1e-16)),
//         max_iter(params.GetOptionalParam<unsigned int>("MaxIter", 1000)),
//         print_level(params.GetOptionalParam<int>("PrintLevel", -1)) {

//     amg.SetPrintLevel(print_level);
//     SetTol(tol);
//     SetMaxIter(max_iter);
//     SetPrintLevel(print_level);
//     SetPreconditioner(amg);
//   }
//   mfem::HypreBoomerAMG amg;
//   double tol;
//   int max_iter;
//   int print_level;
// };

} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

std::string GetTimeDerivativeName(const std::string &name);

std::vector<std::string>
GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain EM formulation.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator {
public:
  TimeDomainProblemOperator(
      mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
      hephaestus::GridFunctions &gridfunctions, hephaestus::BCMap &bc_map,
      hephaestus::Coefficients &coefficients, hephaestus::Sources &sources,
      mfem::Solver &jacobian_solver, mfem::NewtonSolver &nonlinear_solver)
      : myid_(0), num_procs_(1), pmesh_(&pmesh), _fespaces(fespaces),
        _gridfunctions(gridfunctions), _bc_map(bc_map), _sources(sources),
        _coefficients(coefficients), _jacobian_solver(&jacobian_solver),
        _nonlinear_solver(nonlinear_solver) {
    MPI_Comm_size(pmesh.GetComm(), &num_procs_);
    MPI_Comm_rank(pmesh.GetComm(), &myid_);
  };

  ~TimeDomainProblemOperator(){};

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector &X);
  virtual void ImplicitSolve(const double dt, const mfem::Vector &X,
                             mfem::Vector &dX_dt) override;

  virtual void buildEquationSystemOperator(double dt);

  void
  SetEquationSystem(hephaestus::TimeDependentEquationSystem *equation_system);

  virtual mfem::Solver *getJacobianSolver() { return _jacobian_solver; };
  virtual mfem::NewtonSolver *getNonlinearSolver() {
    return &_nonlinear_solver;
  };

  mfem::Array<int> true_offsets, block_trueOffsets;

  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> trial_var_names;

  std::vector<mfem::ParGridFunction *> trial_variable_time_derivatives,
      trial_variables;

  hephaestus::TimeDependentEquationSystem *_equation_system;

  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  hephaestus::FESpaces &_fespaces;
  hephaestus::GridFunctions &_gridfunctions;
  hephaestus::BCMap &_bc_map;
  hephaestus::Sources &_sources;
  hephaestus::Coefficients &_coefficients;
  mfem::BlockVector trueX, trueRhs;
  mfem::OperatorHandle _equation_system_operator;

private:
  mfem::Solver *_jacobian_solver;
  mfem::NewtonSolver &_nonlinear_solver;
};

} // namespace hephaestus

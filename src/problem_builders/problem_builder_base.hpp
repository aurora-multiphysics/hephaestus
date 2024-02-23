#pragma once
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "gridfunctions.hpp"
#include "inputs.hpp"
#include "sources.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus
{

class Problem
{
public:
  Problem() = default;
  virtual ~Problem();

  std::shared_ptr<mfem::ParMesh> _pmesh{nullptr};
  hephaestus::BCMap _bc_map;
  hephaestus::Coefficients _coefficients;
  hephaestus::AuxSolvers _preprocessors;
  hephaestus::AuxSolvers _postprocessors;
  hephaestus::Sources _sources;
  hephaestus::Outputs _outputs;
  hephaestus::InputParameters _solver_options;

  std::unique_ptr<mfem::ODESolver> _ode_solver{nullptr};
  std::unique_ptr<mfem::BlockVector> _f{nullptr};

  std::shared_ptr<mfem::Solver> _jacobian_preconditioner{nullptr};
  std::shared_ptr<mfem::Solver> _jacobian_solver{nullptr};
  std::shared_ptr<mfem::NewtonSolver> _nonlinear_solver{nullptr};

  hephaestus::FECollections _fecs;
  hephaestus::FESpaces _fespaces;
  hephaestus::GridFunctions _gridfunctions;

  MPI_Comm _comm;
  int _myid;
  int _num_procs;

  [[nodiscard]] virtual hephaestus::EquationSystem * GetEquationSystem() const = 0;
  [[nodiscard]] virtual mfem::Operator * GetOperator() const = 0;
};

class ProblemBuilder
{
private:
  virtual hephaestus::Problem * GetProblem() = 0;

public:
  ProblemBuilder() = default;

  // Virtual destructor required to prevent leaks.
  virtual ~ProblemBuilder() = default;

  void SetMesh(std::shared_ptr<mfem::ParMesh> pmesh);
  void SetFESpaces(hephaestus::FESpaces & fespaces);
  void SetGridFunctions(hephaestus::GridFunctions & gridfunctions);
  void SetBoundaryConditions(hephaestus::BCMap & bc_map);
  void SetAuxSolvers(hephaestus::AuxSolvers & preprocessors);
  void SetPostprocessors(hephaestus::AuxSolvers & postprocessors);
  void SetSources(hephaestus::Sources & sources);
  void SetOutputs(hephaestus::Outputs & outputs);
  void SetSolverOptions(hephaestus::InputParameters & solver_options);
  void SetJacobianPreconditioner(std::shared_ptr<mfem::Solver> preconditioner);
  void SetJacobianSolver(std::shared_ptr<mfem::Solver> solver);
  void SetCoefficients(hephaestus::Coefficients & coefficients);

  void AddFESpace(std::string fespace_name,
                  std::string fec_name,
                  int vdim = 1,
                  int ordering = mfem::Ordering::byNODES);
  void AddGridFunction(std::string gridfunction_name, std::string fespace_name);

  template <class T>
  void AddKernel(std::string var_name, std::shared_ptr<hephaestus::Kernel<T>> kernel)
  {
    GetProblem()->GetEquationSystem()->AddTrialVariableNameIfMissing(var_name);
    GetProblem()->GetEquationSystem()->AddKernel(var_name, std::move(kernel));
  }

  void AddBoundaryCondition(std::string bc_name, std::shared_ptr<hephaestus::BoundaryCondition> bc);
  void AddAuxSolver(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux);
  void AddPostprocessor(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux);
  void AddSource(std::string source_name, std::shared_ptr<hephaestus::Source> source);

  virtual void RegisterFESpaces() = 0;
  virtual void RegisterGridFunctions() = 0;
  virtual void RegisterAuxSolvers() = 0;
  virtual void RegisterCoefficients() = 0;

  virtual void InitializeKernels() = 0;
  virtual void ConstructEquationSystem() = 0;
  virtual void SetOperatorGridFunctions() = 0;
  virtual void ConstructJacobianPreconditioner();
  virtual void ConstructJacobianSolver();
  virtual void ConstructNonlinearSolver();
  virtual void ConstructOperator() = 0;
  virtual void ConstructState() = 0;
  virtual void ConstructTimestepper() = 0;

  void InitializeAuxSolvers();
  void InitializeOutputs();

  /**
   * Call to setup a problem. Similar to "ConstructEquationSystemProblem" in the removed
   * ProblemBuilderSequencer.
   */
  void FinalizeProblem();

protected:
  /// Supported Jacobian solver types.
  enum class SolverType
  {
    HYPRE_PCG,
    HYPRE_GMRES,
    HYPRE_FGMRES,
    HYPRE_AMG,
    SUPER_LU
  };

  /// Structure containing default parameters which can be passed to "ConstructJacobianSolverWithOptions".
  /// These will be used if the user has not supplied their own values.
  struct SolverParams
  {
    double _tolerance;
    double _abs_tolerance;

    unsigned int _max_iteration;

    int _print_level;
    int _k_dim;
  };

  /// Called in "ConstructJacobianSolver". This will create a solver of the chosen type and use the user's input
  /// parameters if they have been provided.
  void ConstructJacobianSolverWithOptions(SolverType type,
                                          SolverParams default_params = {
                                              ._tolerance = 1e-16,
                                              ._abs_tolerance = 1e-16,
                                              ._max_iteration = 1000,
                                              ._print_level = GetGlobalPrintLevel(),
                                              ._k_dim = 10});
};
} // namespace hephaestus

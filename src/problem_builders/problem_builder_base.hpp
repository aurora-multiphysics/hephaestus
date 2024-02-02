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

class ProblemSolvers
{

public:
  ProblemSolvers() = default;
  ProblemSolvers(hephaestus::InputParameters solver_options, MPI_Comm comm)
    : _solver_options{std::move(solver_options)}, _comm{std::move(comm)}, _comm_set{true}
  {
  }
  ~ProblemSolvers() = default;

  void
  SetSolver(std::shared_ptr<mfem::Solver> & solver, std::string solve_type, mfem::Operator & A);
  void SetLinearPreconditioner(mfem::Operator & A);
  void SetLinearSolver(mfem::Operator & A);

  void SetComm(MPI_Comm comm);
  void SetSolverOptions(hephaestus::InputParameters solver_options);

  std::shared_ptr<mfem::Solver> _linear_solver{nullptr};
  std::shared_ptr<mfem::Solver> _linear_preconditioner{nullptr};

private:
  MPI_Comm _comm;
  bool _comm_set{false};
  hephaestus::InputParameters _solver_options;
};

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

  ProblemSolvers _solvers;

  hephaestus::FECollections _fecs;
  hephaestus::FESpaces _fespaces;
  hephaestus::GridFunctions _gridfunctions;

  int _myid;
  int _num_procs;

  virtual hephaestus::EquationSystem * GetEquationSystem() = 0;
  virtual mfem::Operator * GetOperator() = 0;
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
  void SetCoefficients(hephaestus::Coefficients & coefficients);

  void AddFESpace(std::string fespace_name,
                  std::string fec_name,
                  int vdim = 1,
                  int ordering = mfem::Ordering::byNODES);
  void AddGridFunction(std::string gridfunction_name, std::string fespace_name);

  template <class T>
  void AddKernel(std::string var_name, std::shared_ptr<hephaestus::Kernel<T>> kernel)
  {
    GetProblem()->GetEquationSystem()->AddVariableNameIfMissing(var_name);
    GetProblem()->GetEquationSystem()->AddKernel(var_name, std::move(kernel));
  }

  template <class T>
  void AddKernel(std::string var_name, hephaestus::Kernel<T> * kernel)
  {
    GetProblem()->GetEquationSystem()->AddVariableNameIfMissing(var_name);
    GetProblem()->GetEquationSystem()->AddKernel(var_name, std::shared_ptr<Kernel<T>>(kernel));
  }

  void AddBoundaryCondition(std::string bc_name, std::shared_ptr<hephaestus::BoundaryCondition> bc);
  void AddAuxSolver(std::string auxsolver_name, hephaestus::AuxSolver * aux, bool own_data);
  void AddAuxSolver(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux);
  void AddPostprocessor(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux);
  void AddSource(std::string source_name, std::shared_ptr<hephaestus::Source> source);

  virtual void RegisterFESpaces() = 0;
  virtual void RegisterGridFunctions() = 0;
  virtual void RegisterAuxSolvers() = 0;
  virtual void RegisterCoefficients() = 0;

  virtual void InitializeKernels() = 0;
  virtual void ConstructEquationSystem() = 0;
  virtual void ConstructOperator() = 0;
  virtual void ConstructState() = 0;
  virtual void ConstructSolver() = 0;

  void InitializeAuxSolvers();
  void InitializeOutputs();
};

class ProblemBuildSequencer
{
  /**
   * @var Builder
   */
private:
  hephaestus::ProblemBuilder * _problem_builder{nullptr};

  /**
   * The ProblemBuildSequencer works with any builder instance that the client
   * code passes to it. This way, the client code may alter the final type of
   * the newly assembled product.
   */

public:
  ProblemBuildSequencer(hephaestus::ProblemBuilder * problem_builder)
    : _problem_builder{problem_builder}
  {
  }

  /**
   * The ProblemBuildSequencer can construct variations of Problems using the
   * same building steps.
   */
  void ConstructOperatorProblem()
  {
    // SteadyStateProblem
    _problem_builder->RegisterFESpaces();
    _problem_builder->RegisterGridFunctions();
    _problem_builder->RegisterAuxSolvers();
    _problem_builder->RegisterCoefficients();
    _problem_builder->InitializeKernels();
    _problem_builder->ConstructOperator();
    _problem_builder->ConstructState();
    _problem_builder->InitializeAuxSolvers();
    _problem_builder->InitializeOutputs();
  }
  void ConstructEquationSystemProblem()
  {
    _problem_builder->RegisterFESpaces();
    _problem_builder->RegisterGridFunctions();
    _problem_builder->RegisterAuxSolvers();
    _problem_builder->RegisterCoefficients();
    _problem_builder->ConstructEquationSystem();
    _problem_builder->InitializeKernels();
    _problem_builder->ConstructOperator();
    _problem_builder->ConstructState();
    _problem_builder->ConstructSolver();
    _problem_builder->InitializeAuxSolvers();
    _problem_builder->InitializeOutputs();
  }
};

} // namespace hephaestus

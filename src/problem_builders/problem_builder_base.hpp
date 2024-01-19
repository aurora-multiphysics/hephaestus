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

  std::shared_ptr<mfem::ParMesh> pmesh{nullptr};
  hephaestus::BCMap bc_map;
  hephaestus::Coefficients coefficients;
  hephaestus::AuxSolvers preprocessors;
  hephaestus::AuxSolvers postprocessors;
  hephaestus::Sources sources;
  hephaestus::Outputs outputs;
  hephaestus::InputParameters solver_options;

  std::unique_ptr<mfem::ODESolver> ode_solver{nullptr};
  std::unique_ptr<mfem::BlockVector> F{nullptr};

  hephaestus::FECollections fecs;
  hephaestus::FESpaces fespaces;
  hephaestus::GridFunctions gridfunctions;

  int myid_;
  int num_procs_;

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
  void AddKernel(std::string var_name, std::unique_ptr<hephaestus::Kernel<T>> && kernel)
  {
    GetProblem()->GetEquationSystem()->addVariableNameIfMissing(var_name);
    GetProblem()->GetEquationSystem()->addKernel(var_name, kernel);
  }

  template <class T>
  void AddKernel(std::string var_name, hephaestus::Kernel<T> * kernel)
  {
    GetProblem()->GetEquationSystem()->addVariableNameIfMissing(var_name);
    GetProblem()->GetEquationSystem()->addKernel(var_name, std::unique_ptr<Kernel<T>>(kernel));
  }

  void AddBoundaryCondition(std::string bc_name, hephaestus::BoundaryCondition * bc, bool own_data);
  void AddAuxSolver(std::string auxsolver_name, hephaestus::AuxSolver * aux, bool own_data);
  void AddPostprocessor(std::string auxsolver_name, hephaestus::AuxSolver * aux, bool own_data);
  void AddSource(std::string source_name, hephaestus::Source * source, bool own_data);

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

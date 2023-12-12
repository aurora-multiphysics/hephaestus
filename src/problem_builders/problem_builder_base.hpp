#pragma once
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "gridfunctions.hpp"
#include "inputs.hpp"
#include "sources.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

class Problem {
public:
  std::shared_ptr<mfem::ParMesh> pmesh;
  hephaestus::BCMap bc_map;
  hephaestus::Coefficients coefficients;
  hephaestus::AuxSolvers preprocessors;
  hephaestus::AuxSolvers postprocessors;
  hephaestus::Sources sources;
  hephaestus::Outputs outputs;

  mfem::ODESolver *ode_solver;
  mfem::BlockVector *F;

  std::shared_ptr<mfem::Solver> jacobian_preconditioner;
  std::shared_ptr<mfem::Solver> jacobian_solver;
  std::shared_ptr<mfem::NewtonSolver> nonlinear_solver;

  hephaestus::FECollections fecs;
  hephaestus::FESpaces fespaces;
  hephaestus::GridFunctions gridfunctions;
  MPI_Comm comm;
  int myid_;
  int num_procs_;

  Problem() = default;

  virtual hephaestus::EquationSystem *GetEquationSystem() = 0;
  virtual mfem::Operator *GetOperator() = 0;
};

class ProblemBuilder {
private:
  virtual hephaestus::Problem *GetProblem() = 0;

public:
  ProblemBuilder(){};

  void SetMesh(std::shared_ptr<mfem::ParMesh> pmesh);
  void SetFESpaces(hephaestus::FESpaces &fespaces);
  void SetGridFunctions(hephaestus::GridFunctions &gridfunctions);
  void SetBoundaryConditions(hephaestus::BCMap &bc_map);
  void SetAuxSolvers(hephaestus::AuxSolvers &preprocessors);
  void SetPostprocessors(hephaestus::AuxSolvers &postprocessors);
  void SetSources(hephaestus::Sources &sources);
  void SetOutputs(hephaestus::Outputs &outputs);
  void SetJacobianPreconditioner(std::shared_ptr<mfem::Solver> preconditioner);
  void SetJacobianSolver(std::shared_ptr<mfem::Solver> solver);
  void SetCoefficients(hephaestus::Coefficients &coefficients);

  void AddFESpace(std::string fespace_name, std::string fec_name, int vdim = 1,
                  int ordering = mfem::Ordering::byNODES);
  void AddGridFunction(std::string gridfunction_name, std::string fespace_name);
  template <class T>
  void AddKernel(std::string var_name, hephaestus::Kernel<T> *kernel) {
    GetProblem()->GetEquationSystem()->addTrialVariableNameIfMissing(var_name);
    GetProblem()->GetEquationSystem()->addKernel(var_name, kernel);
  };
  void AddBoundaryCondition(std::string bc_name,
                            hephaestus::BoundaryCondition *bc, bool own_data);
  void AddAuxSolver(std::string auxsolver_name, hephaestus::AuxSolver *aux,
                    bool own_data);
  void AddPostprocessor(std::string auxsolver_name, hephaestus::AuxSolver *aux,
                        bool own_data);
  void AddSource(std::string source_name, hephaestus::Source *source,
                 bool own_data);

  virtual void RegisterFESpaces() = 0;
  virtual void RegisterGridFunctions() = 0;
  virtual void RegisterAuxSolvers() = 0;
  virtual void RegisterCoefficients() = 0;

  virtual void InitializeKernels() = 0;
  virtual void ConstructEquationSystem() = 0;
  virtual void ConstructJacobianPreconditioner();
  virtual void ConstructJacobianSolver();
  virtual void ConstructNonlinearSolver();
  virtual void ConstructOperator() = 0;
  virtual void ConstructState() = 0;
  virtual void ConstructTimestepper() = 0;

  void InitializeAuxSolvers();
  void InitializeOutputs();
};

class ProblemBuildSequencer {
  /**
   * @var Builder
   */
private:
  std::unique_ptr<hephaestus::ProblemBuilder> problem_builder;
  /**
   * The ProblemBuildSequencer works with any builder instance that the client
   * code passes to it. This way, the client code may alter the final type of
   * the newly assembled product.
   */

public:
  ProblemBuildSequencer(hephaestus::ProblemBuilder *problem_builder_)
      : problem_builder(problem_builder_) {}

  /**
   * The ProblemBuildSequencer can construct variations of Problems using the
   * same building steps.
   */
  void ConstructOperatorProblem() {
    // SteadyStateProblem
    problem_builder->RegisterFESpaces();
    problem_builder->RegisterGridFunctions();
    problem_builder->RegisterAuxSolvers();
    problem_builder->RegisterCoefficients();
    problem_builder->InitializeKernels();
    problem_builder->ConstructJacobianPreconditioner();
    problem_builder->ConstructJacobianSolver();
    problem_builder->ConstructNonlinearSolver();
    problem_builder->ConstructOperator();
    problem_builder->ConstructState();
    problem_builder->InitializeAuxSolvers();
    problem_builder->InitializeOutputs();
  }
  void ConstructEquationSystemProblem() {
    problem_builder->RegisterFESpaces();
    problem_builder->RegisterGridFunctions();
    problem_builder->RegisterAuxSolvers();
    problem_builder->RegisterCoefficients();
    problem_builder->ConstructEquationSystem();
    problem_builder->InitializeKernels();
    problem_builder->ConstructJacobianPreconditioner();
    problem_builder->ConstructJacobianSolver();
    problem_builder->ConstructNonlinearSolver();
    problem_builder->ConstructOperator();
    problem_builder->ConstructState();
    problem_builder->ConstructTimestepper();
    problem_builder->InitializeAuxSolvers();
    problem_builder->InitializeOutputs();
  }
};

} // namespace hephaestus

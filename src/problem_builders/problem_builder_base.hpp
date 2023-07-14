#pragma once
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "gridfunctions.hpp"
#include "inputs.hpp"
#include "postprocessors.hpp"
#include "sources.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

/** Class that store objects that define the set of PDEs to be solved. */
class Problem {
public:
  std::shared_ptr<mfem::ParMesh> pmesh;  /**< Shared pointer to mesh */
  hephaestus::BCMap bc_map;              /**< Boundary conditions to apply */
  hephaestus::Coefficients coefficients; /**< Global store of coefficients */
  hephaestus::AuxSolvers preprocessors; /**< AuxSolvers to apply before solve */
  hephaestus::AuxSolvers postprocessors; /**< AuxSolvers to apply after solve */
  hephaestus::Sources sources;           /**< Source terms to apply to PDE */
  hephaestus::Outputs outputs; /**< Stores objects controlling outputs */
  hephaestus::InputParameters solver_options; /**< Options for PDE solver(s) */

  mfem::ODESolver *ode_solver; /**< ODE solver used for time evolution */
  mfem::BlockVector *F;        /**< State vector storing DOFs of solution */

  hephaestus::FECollections fecs;          /**< Finite element collections */
  hephaestus::FESpaces fespaces;           /**< Finite element spaces */
  hephaestus::GridFunctions gridfunctions; /**< Finite element variables */
  int myid_;                               /**< ID of the current MPI rank */
  int num_procs_;                          /**< Total number of MPI tasks */

  Problem() = default; /**< Default constructor */

  /** Pure virtual method to return a pointer to the EquationSystem defining the
   * PDE
   */
  virtual hephaestus::EquationSystem *GetEquationSystem() = 0;

  /** Pure virtual method to return a pointer to the mfem::Operator used to
   * solve the PDE */
  virtual mfem::Operator *GetOperator() = 0;
};

/** Builder class to set up an initialised hephaestus::Problem object. */
class ProblemBuilder {
private:
  /** Pure virtual method to return a pointer to the currently stored Problem.
   */
  virtual hephaestus::Problem *GetProblem() = 0;

public:
  ProblemBuilder() = default; /**< Default constructor */

  /** Set the Mesh the PDE will be solved on. */
  void SetMesh(std::shared_ptr<mfem::ParMesh> pmesh);

  /** Set the FESpaces used by GridFunctions in the Problem. */
  void SetFESpaces(hephaestus::FESpaces &fespaces);

  /** Set the FE variables used in the Problem. Should refer to FESpaces that
   * have previously been set by SetFESpaces. */
  void SetGridFunctions(hephaestus::GridFunctions &gridfunctions);

  /** Set BoundaryConditions applied in the Problem. */
  void SetBoundaryConditions(hephaestus::BCMap &bc_map);

  /** Set AuxSolvers to be evaluated before the main solve. */
  void SetAuxSolvers(hephaestus::AuxSolvers &preprocessors);

  /** Set AuxSolvers to be evaluated after the main solve. */
  void SetPostprocessors(hephaestus::AuxSolvers &postprocessors);

  /** Set AuxSolvers to be evaluated after the main solve. */
  void SetSources(hephaestus::Sources &sources);

  /** Set Outputs defining where solution data should be exported to. */
  void SetOutputs(hephaestus::Outputs &outputs);

  /** Set options for the PDE solver. */
  void SetSolverOptions(hephaestus::InputParameters &solver_options);

  /** Set coefficients used in Problem solve. */
  void SetCoefficients(hephaestus::Coefficients &coefficients);

  /** Add one mfem::ParFiniteElementSpace to the Problem, with associated
   * mfem::FiniteElementCollection.
   * @param fespace_name Name of new FESpace to add. Should be unique.
   * @param fec_name Name of mfem::FiniteElementCollection passed to
   * mfem::FiniteElementCollection::New to add new FECollection.
   * @param vdim Dimension of the coefficients of shape functions used.
   * @param ordering Enum controlling ordering of DOFs in FESpace
   */
  void AddFESpace(std::string fespace_name, std::string fec_name, int vdim = 1,
                  int ordering = mfem::Ordering::byNODES);

  /** Add one mfem::GridFunction to the Problem.
   * @param gridfunction_name Name of new GridFunction to add. Should be unique.
   * @param fespace_name Name of FESpace stored in Problem->FESpaces that the
   * new GridFunction is defined on.
   */
  void AddGridFunction(std::string gridfunction_name, std::string fespace_name);

  /** Add one hephaestus::Kernel to the Problem.
   * @param var_name Name of GridFunction corresponding to the test variable
   * this kernel will be applied to.
   * @param kernel Pointer to the kernel to add to the Problem.
   */
  template <class T>
  void AddKernel(std::string var_name, hephaestus::Kernel<T> *kernel) {
    this->GetProblem()->GetEquationSystem()->addVariableNameIfMissing(var_name);
    this->GetProblem()->GetEquationSystem()->addKernel(var_name, kernel);
  };

  /** Pure virtual function to register default FESpaces for this Problem. */
  virtual void RegisterFESpaces() = 0;

  /** Pure virtual function to register default GridFunctions for this Problem.
   */
  virtual void RegisterGridFunctions() = 0;

  /** Pure virtual function to register default AuxSolvers for this Problem. */
  virtual void RegisterAuxSolvers() = 0;

  /** Pure virtual function to register required Coefficients for this Problem.
   */
  virtual void RegisterCoefficients() = 0;

  /** Pure virtual function to intialize all stored kernels (including sources)
   * on the EquationSystem associated with this problem. */
  virtual void InitializeKernels() = 0;

  /** Pure virtual function to construct EquationSystem storing all linear and
   * bilinear forms associated with this Problem. */
  virtual void ConstructEquationSystem() = 0;

  /** Pure virtual function to construct mfem::Operator used to solve Problem.
   */
  virtual void ConstructOperator() = 0;

  /** Pure virtual function to construct the initial state vector for the
   * problem. */
  virtual void ConstructState() = 0;

  /** Pure virtual function to construct the ODE solver for time evolution. */
  virtual void ConstructSolver() = 0;

  /** Initialize Postprocessors for this problem. */
  void InitializePostprocessors();
};

/** Director class for ProblemBuilder controlling order of build operations. */
class ProblemBuildSequencer {
private:
  /** Unique pointer to problem_builder managed by this sequencer. */
  std::unique_ptr<hephaestus::ProblemBuilder> problem_builder;

public:
  /** Unique pointer to problem_builder controlled by this sequencer.
   * @param problem_builder Pointer to ProblemBuilder to be managed by this
   * sequencer.
   */
  ProblemBuildSequencer(hephaestus::ProblemBuilder *problem_builder_)
      : problem_builder(problem_builder_) {}

  /**
   * Construct a Problem defined by just an Operator to use in the linear
   * solve. Appropriate for steady state or frequency domain solves.
   */
  void ConstructOperatorProblem() {
    this->problem_builder->RegisterFESpaces();
    this->problem_builder->RegisterGridFunctions();
    this->problem_builder->RegisterAuxSolvers();
    this->problem_builder->RegisterCoefficients();
    this->problem_builder->InitializeKernels();
    this->problem_builder->ConstructOperator();
    this->problem_builder->ConstructState();
    this->problem_builder->InitializePostprocessors();
  }

  /**
   * Construct a Problem associated with an EquationSystem.
   */
  void ConstructEquationSystemProblem() {
    this->problem_builder->RegisterFESpaces();
    this->problem_builder->RegisterGridFunctions();
    this->problem_builder->RegisterAuxSolvers();
    this->problem_builder->RegisterCoefficients();
    this->problem_builder->ConstructEquationSystem();
    this->problem_builder->InitializeKernels();
    this->problem_builder->ConstructOperator();
    this->problem_builder->ConstructState();
    this->problem_builder->ConstructSolver();
    this->problem_builder->InitializePostprocessors();
  }
};

} // namespace hephaestus

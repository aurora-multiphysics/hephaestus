#pragma once
#include "auxsolvers.hpp"
#include "inputs.hpp"
#include "postprocessors.hpp"
#include "sources.hpp"
#include "variables.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

class Problem {
public:
  std::shared_ptr<mfem::ParMesh> pmesh;
  hephaestus::BCMap bc_map;
  hephaestus::DomainProperties domain_properties;
  hephaestus::AuxSolvers preprocessors;
  hephaestus::AuxSolvers postprocessors;
  hephaestus::Sources sources;
  hephaestus::Outputs outputs;
  std::map<std::string, mfem::DataCollection *> data_collections;
  hephaestus::InputParameters solver_options;

  mfem::ODESolver *ode_solver;
  mfem::BlockVector *F;

  mfem::NamedFieldsMap<mfem::FiniteElementCollection> fecs;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> gridfunctions;
  int myid_;
  int num_procs_;

  Problem() = default;
  explicit Problem(const hephaestus::InputParameters &params);

  virtual hephaestus::Formulation *GetFormulation() = 0;
  virtual hephaestus::EquationSystem *GetEquationSystem() = 0;
  virtual mfem::Operator *GetOperator() = 0;
};

class ProblemBuilder {
private:
  virtual hephaestus::Problem *GetProblem() = 0;

public:
  ProblemBuilder(){};
  ProblemBuilder(const hephaestus::InputParameters &params){};

  virtual void SetMesh(std::shared_ptr<mfem::ParMesh> pmesh) {
    this->GetProblem()->pmesh = pmesh;
  };

  virtual void SetFESpaces(hephaestus::FESpaces &fespaces) {
    this->GetProblem()->fespaces = fespaces;
  };

  virtual void SetGridFunctions(hephaestus::GridFunctions &gridfunctions) {
    this->GetProblem()->gridfunctions = gridfunctions;
  };

  virtual void SetBoundaryConditions(hephaestus::BCMap &bc_map) {
    this->GetProblem()->bc_map = bc_map;
  };

  virtual void SetAuxSolvers(hephaestus::AuxSolvers &preprocessors) {
    this->GetProblem()->preprocessors = preprocessors;
  };

  virtual void SetPostprocessors(hephaestus::AuxSolvers &postprocessors) {
    this->GetProblem()->postprocessors = postprocessors;
  };

  virtual void SetSources(hephaestus::Sources &sources) {
    this->GetProblem()->sources = sources;
  };

  virtual void SetOutputs(hephaestus::Outputs &outputs) {
    this->GetProblem()->outputs = outputs;
  };

  virtual void SetSolverOptions(hephaestus::InputParameters &solver_options) {
    this->GetProblem()->solver_options = solver_options;
  };

  virtual void
  SetCoefficients(hephaestus::DomainProperties &domain_properties) {
    this->GetProblem()->domain_properties = domain_properties;
  };

  void AddFESpace(std::string fespace_name, std::string fec_name, int vdim = 1,
                  int ordering = mfem::Ordering::byNODES) {
    mfem::FiniteElementCollection *fec =
        mfem::FiniteElementCollection::New(fec_name.c_str());
    this->GetProblem()->fecs.Register(fec_name, fec, true);

    mfem::ParMesh *pmesh = this->GetProblem()->pmesh.get();
    if (pmesh == NULL) {
      MFEM_ABORT("ParMesh not found when trying to add " << fespace_name
                                                         << " to fespaces.");
    }

    mfem::ParFiniteElementSpace *pfes = new mfem::ParFiniteElementSpace(
        this->GetProblem()->pmesh.get(), fec, vdim, ordering);

    this->GetProblem()->fespaces.Register(fespace_name, pfes, true);
  };

  void AddGridFunction(std::string gridfunction_name,
                       std::string fespace_name) {
    mfem::ParFiniteElementSpace *fespace(
        this->GetProblem()->fespaces.Get(fespace_name));
    if (fespace == NULL) {
      MFEM_ABORT("FESpace "
                 << fespace_name << " not found in fespaces when trying to add "
                 << gridfunction_name
                 << " associated with it into gridfunctions. Please add "
                 << fespace_name
                 << " to fespaces before adding this gridfunction.");
    }
    mfem::ParGridFunction *gridfunc = new mfem::ParGridFunction(fespace);
    *gridfunc = 0.0;

    this->GetProblem()->gridfunctions.Register(gridfunction_name, gridfunc,
                                               true);
  };

  template <typename T>
  void AddKernel(std::string var_name, hephaestus::Kernel<T> *kernel) {
    this->GetProblem()->GetEquationSystem()->addVariableNameIfMissing(var_name);
    this->GetProblem()->GetEquationSystem()->addKernel(var_name, kernel);
  };

  virtual void RegisterFESpaces() = 0;

  virtual void RegisterGridFunctions() = 0;

  virtual void RegisterAuxSolvers() = 0;

  virtual void RegisterCoefficients() = 0;

  virtual void InitializeKernels() = 0;

  virtual void ConstructEquationSystem() = 0;

  virtual void ConstructOperator() = 0;

  virtual void ConstructState() = 0;

  virtual void ConstructSolver() = 0;

  virtual void InitializePostprocessors() = 0;
}; // namespace hephaestus

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
    this->problem_builder->RegisterFESpaces();
    this->problem_builder->RegisterGridFunctions();
    this->problem_builder->RegisterAuxSolvers();
    this->problem_builder->RegisterCoefficients();
    this->problem_builder->InitializeKernels();
    this->problem_builder->ConstructOperator();
    this->problem_builder->ConstructState();
    this->problem_builder->InitializePostprocessors();
  }
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

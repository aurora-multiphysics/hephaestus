#pragma once
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "hephaestus_containers.hpp"
#include "inputs.hpp"
#include "sources.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus
{

/// Forwards declaration.
class ProblemOperatorBase;

/// Base problem class.
class Problem
{
public:
  Problem() = default;
  virtual ~Problem() = default;

  std::shared_ptr<mfem::ParMesh> _pmesh{nullptr};
  hephaestus::BCMap _bc_map;
  hephaestus::Coefficients _coefficients;
  hephaestus::AuxSolvers _preprocessors;
  hephaestus::AuxSolvers _postprocessors;
  hephaestus::Sources _sources;
  hephaestus::Outputs _outputs;

  hephaestus::FECollections _fecs;
  hephaestus::FESpaces _fespaces;
  hephaestus::GridFunctions _gridfunctions;

  MPI_Comm _comm;
  int _myid;
  int _num_procs;

  /// Returns a pointer to the operator. See derived classes.
  [[nodiscard]] virtual hephaestus::ProblemOperatorBase * GetOperator() const = 0;

  /// Call to update on mesh change.
  virtual void Update();
};

/// ProblemBuilder base class.
class ProblemBuilder
{
public:
  /// NB: delete empty constructor to allow only derived classes to be constructed.
  ProblemBuilder() = delete;

  // Virtual destructor required to prevent leaks.
  virtual ~ProblemBuilder()
  {
    if (_problem != nullptr)
    {
      delete _problem;
    }
  }

  void SetMesh(std::shared_ptr<mfem::ParMesh> pmesh);
  void SetFESpaces(hephaestus::FESpaces & fespaces);
  void SetGridFunctions(hephaestus::GridFunctions & gridfunctions);
  void SetBoundaryConditions(hephaestus::BCMap & bc_map);
  void SetAuxSolvers(hephaestus::AuxSolvers & preprocessors);
  void SetPostprocessors(hephaestus::AuxSolvers & postprocessors);
  void SetSources(hephaestus::Sources & sources);
  void SetOutputs(hephaestus::Outputs & outputs);
  void SetCoefficients(hephaestus::Coefficients & coefficients);

  void AddFESpace(std::string fespace_name,
                  std::string fec_name,
                  int vdim = 1,
                  int ordering = mfem::Ordering::byNODES);
  void AddGridFunction(std::string gridfunction_name, std::string fespace_name);

  void AddBoundaryCondition(std::string bc_name, std::shared_ptr<hephaestus::BoundaryCondition> bc);
  void AddAuxSolver(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux);
  void AddPostprocessor(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux);
  void AddSource(std::string source_name, std::shared_ptr<hephaestus::Source> source);

  virtual void RegisterFESpaces() = 0;
  virtual void RegisterGridFunctions() = 0;
  virtual void RegisterAuxSolvers() = 0;
  virtual void RegisterCoefficients() = 0;

  virtual void ConstructOperator() = 0;

  virtual void InitializeSources();

  void InitializeAuxSolvers();
  void InitializeOutputs();

  virtual void InitializeOperator();

  /// @brief Call @a FinalizeProblem to setup a problem.
  /// @param build_operator Skips @a ConstructOperator step if false. Set this to false if the problem
  /// operator has already been constructed earlier to avoid rebuilding it.
  void FinalizeProblem(bool build_operator = true);

protected:
  /// Protected constructor. Derived classes must call this constructor.
  ProblemBuilder(hephaestus::Problem * problem) : _problem{problem} {}

  /// Overridden in derived classes.
  [[nodiscard]] virtual hephaestus::Problem * GetProblem() const = 0;

  /// Helper template getter with safety check.
  template <class TDerivedProblem>
  [[nodiscard]] TDerivedProblem * GetProblem() const
  {
    if (!_problem)
    {
      MFEM_ABORT("hephaestus::Problem instance is NULL.");
    }

    return static_cast<TDerivedProblem *>(_problem);
  }

  /// Helper template for returning a unique pointer of the correct class.
  /// Sets _problem = nullptr to avoid double-free.
  template <class TDerivedProblem>
  [[nodiscard]] std::unique_ptr<TDerivedProblem> ReturnProblem()
  {
    auto * problem = GetProblem<TDerivedProblem>();
    _problem = nullptr; // Avoid double-free.

    return std::unique_ptr<TDerivedProblem>(problem);
  }

  /// Coefficient used in some derived classes.
  mfem::ConstantCoefficient _one_coef{1.0};

private:
  hephaestus::Problem * _problem{nullptr};
};

/// Interface for problem builders that are constructing problems with an equation system.
class EquationSystemProblemBuilderInterface
{
public:
  EquationSystemProblemBuilderInterface() = default;
  virtual ~EquationSystemProblemBuilderInterface() = default;

  /// Add a kernel to the problem's equation system.
  template <class T>
  void AddKernel(std::string var_name, std::shared_ptr<hephaestus::Kernel<T>> kernel)
  {
    GetEquationSystem()->AddTrialVariableNameIfMissing(var_name);
    GetEquationSystem()->AddKernel(var_name, std::move(kernel));
  }

protected:
  /// Implemented in derived classes. Returns a pointer to the problem operator's equation system.
  [[nodiscard]] virtual hephaestus::EquationSystem * GetEquationSystem() const = 0;
};

} // namespace hephaestus
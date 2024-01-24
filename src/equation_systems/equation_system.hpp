#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"
#include "kernel_base.hpp"
#include "named_fields_map.hpp"
#include "sources.hpp"

namespace hephaestus
{

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class EquationSystem
{
public:
  using ParBilinearFormKernel = hephaestus::Kernel<mfem::ParBilinearForm>;
  using ParLinearFormKernel = hephaestus::Kernel<mfem::ParLinearForm>;
  using ParNonlinearFormKernel = hephaestus::Kernel<mfem::ParNonlinearForm>;
  using ParMixedBilinearFormKernel = hephaestus::Kernel<mfem::ParMixedBilinearForm>;

  EquationSystem() = default;
  EquationSystem(const hephaestus::InputParameters & params);

  virtual ~EquationSystem();

  // Names of all gridfunctions corresponding to gridfunctions. This may differ
  // from test_var_names when test gridfunctions include time derivatives.
  std::vector<std::string> _var_names;
  // Names of all test gridfunctions with kernels in this equation system.
  std::vector<std::string> _test_var_names;
  std::vector<mfem::ParFiniteElementSpace *> _test_pfespaces;

  // Components of weak form. // Named according to test variable
  hephaestus::NamedFieldsMap<mfem::ParBilinearForm> _blfs;
  hephaestus::NamedFieldsMap<mfem::ParLinearForm> _lfs;
  hephaestus::NamedFieldsMap<mfem::ParNonlinearForm> _nlfs;
  hephaestus::NamedFieldsMap<hephaestus::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _mblfs; // named according to trial variable

  // add test variable to EquationSystem;
  virtual void AddTestVariableNameIfMissing(const std::string & test_var_name);
  virtual void AddVariableNameIfMissing(const std::string & var_name);

  // Add kernels.
  void AddKernel(const std::string & test_var_name,
                 std::unique_ptr<ParBilinearFormKernel> && blf_kernel);

  void AddKernel(const std::string & test_var_name,
                 std::unique_ptr<ParLinearFormKernel> && lf_kernel);

  void AddKernel(const std::string & test_var_name,
                 std::unique_ptr<ParNonlinearFormKernel> && nlf_kernel);

  void AddKernel(const std::string & trial_var_name,
                 const std::string & test_var_name,
                 std::unique_ptr<ParMixedBilinearFormKernel> && mblf_kernel);

  virtual void ApplyBoundaryConditions(hephaestus::BCMap & bc_map);

  // override to add kernels
  virtual void AddKernels() {}

  // Build forms
  virtual void Init(hephaestus::GridFunctions & gridfunctions,
                    const hephaestus::FESpaces & fespaces,
                    hephaestus::BCMap & bc_map,
                    hephaestus::Coefficients & coefficients);
  virtual void BuildLinearForms(hephaestus::BCMap & bc_map, hephaestus::Sources & sources);
  virtual void BuildBilinearForms();
  virtual void BuildMixedBilinearForms();
  virtual void BuildEquationSystem(hephaestus::BCMap & bc_map, hephaestus::Sources & sources);

  // Form linear system, with essential boundary conditions accounted for
  virtual void FormLinearSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);

  // Update variable from solution vector after solve
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  hephaestus::GridFunctions & gridfunctions);

protected:
  bool VectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;

  // gridfunctions for setting Dirichlet BCs
  std::vector<mfem::Array<int>> _ess_tdof_lists;
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _xs;

  mfem::Array2D<mfem::HypreParMatrix *> _h_blocks;

  // Arrays to store kernels to act on each component of weak form. Named
  // according to test variable
  hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParBilinearFormKernel>>> _blf_kernels_map;

  hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParLinearFormKernel>>> _lf_kernels_map;

  hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParNonlinearFormKernel>>> _nlf_kernels_map;

  hephaestus::NamedFieldsMap<
      hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParMixedBilinearFormKernel>>>>
      _mblf_kernels_map_map;
};

/*
Class to store weak form components for time dependent PDEs
*/
class TimeDependentEquationSystem : public EquationSystem
{
public:
  TimeDependentEquationSystem(const hephaestus::InputParameters & params);
  ~TimeDependentEquationSystem() override = default;

  static std::string GetTimeDerivativeName(std::string name)
  {
    return std::string("d") + name + std::string("_dt");
  }
  void AddVariableNameIfMissing(const std::string & var_name) override;

  virtual void SetTimeStep(double dt);
  virtual void UpdateEquationSystem(hephaestus::BCMap & bc_map, hephaestus::Sources & sources);
  mfem::ConstantCoefficient _dt_coef; // Coefficient for timestep scaling
  std::vector<std::string> _var_time_derivative_names;
};

} // namespace hephaestus

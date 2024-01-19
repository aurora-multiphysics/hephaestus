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
  typedef hephaestus::Kernel<mfem::ParBilinearForm> ParBilinearFormKernel;
  typedef hephaestus::Kernel<mfem::ParLinearForm> ParLinearFormKernel;
  typedef hephaestus::Kernel<mfem::ParNonlinearForm> ParNonlinearFormKernel;
  typedef hephaestus::Kernel<mfem::ParMixedBilinearForm> ParMixedBilinearFormKernel;

  EquationSystem(){};
  EquationSystem(const hephaestus::InputParameters & params);

  virtual ~EquationSystem();

  // Names of all gridfunctions corresponding to gridfunctions. This may differ
  // from test_var_names when test gridfunctions include time derivatives.
  std::vector<std::string> var_names;
  // Names of all test gridfunctions with kernels in this equation system.
  std::vector<std::string> test_var_names;
  std::vector<mfem::ParFiniteElementSpace *> test_pfespaces;

  // Components of weak form. // Named according to test variable
  hephaestus::NamedFieldsMap<mfem::ParBilinearForm> blfs;
  hephaestus::NamedFieldsMap<mfem::ParLinearForm> lfs;
  hephaestus::NamedFieldsMap<mfem::ParNonlinearForm> nlfs;
  hephaestus::NamedFieldsMap<hephaestus::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      mblfs; // named according to trial variable

  // add test variable to EquationSystem;
  virtual void addTestVariableNameIfMissing(const std::string & test_var_name);
  virtual void addVariableNameIfMissing(const std::string & var_name);

  // Add kernels.
  void addKernel(const std::string & test_var_name,
                 std::unique_ptr<ParBilinearFormKernel> && blf_kernel);

  void addKernel(const std::string & test_var_name,
                 std::unique_ptr<ParLinearFormKernel> && lf_kernel);

  void addKernel(const std::string & test_var_name,
                 std::unique_ptr<ParNonlinearFormKernel> && nlf_kernel);

  void addKernel(const std::string & trial_var_name,
                 const std::string & test_var_name,
                 std::unique_ptr<ParMixedBilinearFormKernel> && mblf_kernel);

  virtual void applyBoundaryConditions(hephaestus::BCMap & bc_map);

  // override to add kernels
  virtual void addKernels(){};

  // Build forms
  virtual void Init(hephaestus::GridFunctions & gridfunctions,
                    const hephaestus::FESpaces & fespaces,
                    hephaestus::BCMap & bc_map,
                    hephaestus::Coefficients & coefficients);
  virtual void buildLinearForms(hephaestus::BCMap & bc_map, hephaestus::Sources & sources);
  virtual void buildBilinearForms();
  virtual void buildMixedBilinearForms();
  virtual void buildEquationSystem(hephaestus::BCMap & bc_map, hephaestus::Sources & sources);

  // Form linear system, with essential boundary conditions accounted for
  virtual void FormLinearSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);

  // Update variable from solution vector after solve
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  hephaestus::GridFunctions & gridfunctions);

protected:
  bool vectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;

  // gridfunctions for setting Dirichlet BCs
  std::vector<mfem::Array<int>> ess_tdof_lists;
  std::vector<std::unique_ptr<mfem::ParGridFunction>> xs;

  mfem::Array2D<mfem::HypreParMatrix *> hBlocks;

  // Arrays to store kernels to act on each component of weak form. Named
  // according to test variable
  hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParBilinearFormKernel>>> blf_kernels_map;

  hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParLinearFormKernel>>> lf_kernels_map;

  hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParNonlinearFormKernel>>> nlf_kernels_map;

  hephaestus::NamedFieldsMap<
      hephaestus::NamedFieldsMap<std::vector<std::unique_ptr<ParMixedBilinearFormKernel>>>>
      mblf_kernels_map_map;
};

/*
Class to store weak form components for time dependent PDEs
*/
class TimeDependentEquationSystem : public EquationSystem
{
public:
  TimeDependentEquationSystem(const hephaestus::InputParameters & params);
  ~TimeDependentEquationSystem() override{};

  static std::string GetTimeDerivativeName(std::string name)
  {
    return std::string("d") + name + std::string("_dt");
  }
  virtual void addVariableNameIfMissing(const std::string & var_name) override;

  virtual void setTimeStep(double dt);
  virtual void updateEquationSystem(hephaestus::BCMap & bc_map, hephaestus::Sources & sources);
  mfem::ConstantCoefficient dtCoef; // Coefficient for timestep scaling
  std::vector<std::string> var_time_derivative_names;
};

} // namespace hephaestus

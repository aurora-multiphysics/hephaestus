#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {
/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class EquationSystem {
public:
  EquationSystem(){};
  EquationSystem(const hephaestus::InputParameters &params);

  ~EquationSystem(){};

  // Names of all variables corresponding to variables. This may differ from
  // test_var_names when test variables include time derivatives.
  std::vector<std::string> var_names;
  // Names of all test variables with kernels in this equation system.
  std::vector<std::string> test_var_names;
  std::vector<mfem::ParFiniteElementSpace *> test_pfespaces;

  // Components of weak form. // Named according to test variable
  mfem::NamedFieldsMap<mfem::ParBilinearForm> blfs;
  mfem::NamedFieldsMap<mfem::ParLinearForm> lfs;
  mfem::NamedFieldsMap<mfem::ParNonlinearForm> nlfs;
  mfem::NamedFieldsMap<mfem::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      mblfs; // named according to trial variable

  // add test variable to EquationSystem;
  virtual void addTestVariableNameIfMissing(std::string test_var_name);
  virtual void addVariableNameIfMissing(std::string var_name);

  // Add kernels. EquationSystem takes ownership.
  void addKernel(std::string test_var_name,
                 hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel);
  void addKernel(std::string test_var_name,
                 hephaestus::Kernel<mfem::ParLinearForm> *lf_kernel);
  void addKernel(std::string test_var_name,
                 hephaestus::Kernel<mfem::ParNonlinearForm> *nlf_kernel);
  void addKernel(std::string trial_var_name, std::string test_var_name,
                 hephaestus::Kernel<mfem::ParMixedBilinearForm> *mblf_kernel);
  virtual void applyBoundaryConditions(hephaestus::BCMap &bc_map);

  // override to add kernels
  virtual void addKernels(){};

  // Build forms
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties);
  virtual void buildLinearForms(hephaestus::BCMap &bc_map,
                                hephaestus::Sources &sources);
  virtual void buildBilinearForms();
  virtual void buildMixedBilinearForms();
  virtual void buildEquationSystem(hephaestus::BCMap &bc_map,
                                   hephaestus::Sources &sources);

  // Form linear system, with essential boundary conditions accounted for
  virtual void FormLinearSystem(mfem::OperatorHandle &op,
                                mfem::BlockVector &trueX,
                                mfem::BlockVector &trueRHS);

  // Update variable from solution vector after solve
  virtual void
  RecoverFEMSolution(mfem::BlockVector &trueX,
                     mfem::NamedFieldsMap<mfem::ParGridFunction> &variables);

protected:
  // Variables for setting Dirichlet BCs
  std::vector<mfem::Array<int>> ess_tdof_lists;
  std::vector<mfem::ParGridFunction *> xs;

  mfem::Array2D<mfem::HypreParMatrix *> hBlocks;
  // Arrays to store kernels to act on each component of weak form. Named
  // according to test variable
  mfem::NamedFieldsMap<std::vector<hephaestus::Kernel<mfem::ParBilinearForm> *>>
      blf_kernels_map;
  mfem::NamedFieldsMap<std::vector<hephaestus::Kernel<mfem::ParLinearForm> *>>
      lf_kernels_map;
  mfem::NamedFieldsMap<
      std::vector<hephaestus::Kernel<mfem::ParNonlinearForm> *>>
      nlf_kernels_map;
  mfem::NamedFieldsMap<mfem::NamedFieldsMap<
      std::vector<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>>>
      mblf_kernels_map_map;
};

/*
Class to store weak form components for time dependent PDEs
*/
class TimeDependentEquationSystem : public EquationSystem {
public:
  TimeDependentEquationSystem(const hephaestus::InputParameters &params);
  ~TimeDependentEquationSystem(){};

  static std::string GetTimeDerivativeName(std::string name) {
    return std::string("d") + name + std::string("_dt");
  }
  virtual void addVariableNameIfMissing(std::string var_name) override;

  virtual void setTimeStep(double dt);
  virtual void updateEquationSystem(hephaestus::BCMap &bc_map,
                                    hephaestus::Sources &sources);
  mfem::ConstantCoefficient dtCoef; // Coefficient for timestep scaling
  std::vector<std::string> var_time_derivative_names;
}; // namespace hephaestus

class CurlCurlEquationSystem : public TimeDependentEquationSystem {
public:
  CurlCurlEquationSystem(const hephaestus::InputParameters &params);
  ~CurlCurlEquationSystem(){};
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void addKernels() override;

  std::string h_curl_var_name, alpha_coef_name, beta_coef_name,
      dtalpha_coef_name;
};

class AVEquationSystem : public TimeDependentEquationSystem {
public:
  AVEquationSystem(const hephaestus::InputParameters &params);
  ~AVEquationSystem(){};
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties);
  virtual void addKernels() override;

  std::string a_name, v_name, coupled_variable_name, alpha_coef_name,
      beta_coef_name, dtalpha_coef_name, neg_beta_coef_name;
  mfem::ConstantCoefficient negCoef;
};
} // namespace hephaestus

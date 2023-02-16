#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class WeakForm {
public:
  WeakForm(const hephaestus::InputParameters &params);

  ~WeakForm(){};

  // Name of test variable and pointer to its ParFiniteElementSpace.
  std::string _test_var_name; // TODO: should remove by updating bcmap
  mfem::ParFiniteElementSpace *test_pfes;

  // Components of weak form. All should share the test variable FESpace.
  mfem::ParBilinearForm *blf;
  mfem::ParLinearForm *lf;
  mfem::ParNonlinearForm *nlf;
  mfem::NamedFieldsMap<mfem::ParMixedBilinearForm>
      mblfs; // named according to trial variable

  // Add kernels. WeakForm takes ownership.
  void addKernel(hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel);
  void addKernel(hephaestus::Kernel<mfem::ParLinearForm> *lf_kernel);
  void addKernel(hephaestus::Kernel<mfem::ParNonlinearForm> *nlf_kernel);
  void addKernel(std::string trial_var_name,
                 hephaestus::Kernel<mfem::ParMixedBilinearForm> *mblf_kernel);
  virtual void applyBoundaryConditions(hephaestus::BCMap &bc_map);

  // override to add kernels
  virtual void
  addKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
             const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
             hephaestus::BCMap &bc_map,
             hephaestus::DomainProperties &domain_properties){};

  // Build forms
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties);
  virtual void buildLinearForm(hephaestus::BCMap &bc_map,
                               hephaestus::Sources &sources);
  virtual void buildBilinearForm();
  virtual void buildWeakForm(hephaestus::BCMap &bc_map,
                             hephaestus::Sources &sources);

  // Form linear system, with essential boundary conditions accounted for
  virtual void FormLinearSystem(mfem::HypreParMatrix &A, mfem::Vector &X,
                                mfem::Vector &B);

  // Update variable from solution vector after solve
  virtual void RecoverFEMSolution(mfem::Vector &X,
                                  mfem::ParGridFunction &test_variable);

protected:
  // Variables for setting Dirichlet BCs
  mfem::Array<int> ess_tdof_list;
  mfem::ParGridFunction *x;

  // Arrays to store kernels to act on each component of weak form.
  mfem::Array<hephaestus::Kernel<mfem::ParBilinearForm> *> blf_kernels;
  mfem::Array<hephaestus::Kernel<mfem::ParLinearForm> *> lf_kernels;
  mfem::Array<hephaestus::Kernel<mfem::ParNonlinearForm> *> nlf_kernels;
  mfem::NamedFieldsMap<
      mfem::Array<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>>
      mblf_kernels_map;
};

/*
Class to store weak form components for time dependent PDEs
*/
class TimeDependentWeakForm : public WeakForm {
public:
  TimeDependentWeakForm(const hephaestus::InputParameters &params);
  ~TimeDependentWeakForm(){};
  virtual void setTimeStep(double dt);
  virtual void updateWeakForm(hephaestus::BCMap &bc_map,
                              hephaestus::Sources &sources);
  mfem::ConstantCoefficient dtCoef; // Coefficient for timestep scaling
};

/*
Class to store weak form components for the equation
∇×(α∇×u) + βdu/dt = s0
where
s0 ∈ H(div) source field
u ∈ H(curl)
with corresponding weak form
(α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
- (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
blf(du/dt_{n+1}, u') = lf(u')
blf(u, u') = (βu, u') + (αdt∇×u, ∇×u')
lf(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>
using
u_{n+1} = u_{n} + dt du/dt_{n+1}
*/
class CurlCurlWeakForm : public TimeDependentWeakForm {
public:
  CurlCurlWeakForm(const hephaestus::InputParameters &params);
  ~CurlCurlWeakForm(){};
  virtual void
  addKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
             const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
             hephaestus::BCMap &bc_map,
             hephaestus::DomainProperties &domain_properties) override;
  virtual void setTimeStep(double dt) override;
  virtual void updateWeakForm(hephaestus::BCMap &bc_map,
                              hephaestus::Sources &sources) override;

  std::string coupled_variable_name, alpha_coef_name, beta_coef_name,
      dtalpha_coef_name;
};
} // namespace hephaestus

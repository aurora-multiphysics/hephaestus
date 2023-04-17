// #pragma once
// #include "../common/pfem_extras.hpp"
// #include "inputs.hpp"
// #include "sources.hpp"

// namespace hephaestus {
// /*
// Class to store weak form components (bilinear and linear forms, and
// optionally mixed and nonlinear forms) and build methods
// */
// class ComplexEquationSystem {
// public:
//   ComplexEquationSystem(){};
//   ComplexEquationSystem(const hephaestus::InputParameters &params);

//   virtual ~ComplexEquationSystem();

//   // Names of all variables corresponding to variables. This may differ from
//   // test_var_names when test variables include time derivatives.
//   std::vector<std::string> var_names;
//   // Names of all test variables with kernels in this equation system.
//   std::vector<std::string> test_var_names;
//   std::vector<mfem::ParFiniteElementSpace *> test_pfespaces;

//   // Components of weak form. // Named according to test variable
//   mfem::NamedFieldsMap<mfem::ParBilinearForm> blfs;
//   mfem::NamedFieldsMap<mfem::ParSesquilinearForm> sqlfs;

//   // add test variable to ComplexEquationSystem;
//   virtual void addTestVariableNameIfMissing(std::string test_var_name);
//   virtual void addVariableNameIfMissing(std::string var_name);

//   // Add kernels. ComplexEquationSystem takes ownership.
//   void addKernel(std::string test_var_name,
//                  hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel);
//   void addKernel(std::string test_var_name,
//                  hephaestus::Kernel<mfem::ParSesquilinearForm> *sqlf_kernel);
//   virtual void applyBoundaryConditions(hephaestus::BCMap &bc_map);

//   // override to add kernels
//   virtual void addKernels(){};

//   // Build forms
//   virtual void
//   Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//        const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//        hephaestus::BCMap &bc_map,
//        hephaestus::DomainProperties &domain_properties);
//   virtual void buildBilinearForms(hephaestus::BCMap &bc_map,
//                                   hephaestus::Sources &sources);
//   virtual void buildSesquilinearForms();
//   virtual void buildEquationSystem(hephaestus::BCMap &bc_map,
//                                    hephaestus::Sources &sources);

//   // Form linear system, with essential boundary conditions accounted for
//   virtual void FormLinearSystem(mfem::OperatorHandle &op,
//                                 mfem::BlockVector &trueX,
//                                 mfem::BlockVector &trueRHS);

//   // Update variable from solution vector after solve
//   virtual void
//   RecoverFEMSolution(mfem::BlockVector &trueX,
//                      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables);

// protected:
//   // Variables for setting Dirichlet BCs
//   std::vector<mfem::Array<int>> ess_tdof_lists;
//   std::vector<mfem::ParGridFunction *> xs;

//   mfem::Array2D<mfem::HypreParMatrix *> hBlocks;
//   // Arrays to store kernels to act on each component of weak form. Named
//   // according to test variable
//   mfem::NamedFieldsMap<mfem::Array<hephaestus::Kernel<mfem::ParBilinearForm>
//   *>>
//       blf_kernels_map;
//   mfem::NamedFieldsMap<
//       mfem::Array<hephaestus::Kernel<mfem::ParSesquilinearForm> *>>
//       sqlf_kernels_map;
// };

// } // namespace hephaestus

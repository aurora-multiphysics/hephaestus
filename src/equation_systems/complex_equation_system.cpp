// #include "complex_equation_system.hpp"

// namespace hephaestus {

// ComplexEquationSystem::ComplexEquationSystem(
//     const hephaestus::InputParameters &params)
//     : var_names(), test_var_names(), test_pfespaces(), blfs(), sqlfs(),
//       ess_tdof_lists(), xs() {}

// ComplexEquationSystem::~ComplexEquationSystem() {
//   hBlocks.DeleteAll();
//   blfs.DeleteData(true);
//   sqlfs.DeleteData(true);

//   for (const auto &[test_var_name, blf_kernels] : blf_kernels_map.GetMap()) {
//     blf_kernels->DeleteAll();
//   }
//   blf_kernels_map.DeleteData(true);

//   for (const auto &[test_var_name, sqlf_kernels] : sqlf_kernels_map.GetMap())
//   {
//     sqlf_kernels->DeleteAll();
//   }
//   sqlf_kernels_map.DeleteData(true);
// }

// void ComplexEquationSystem::addVariableNameIfMissing(std::string var_name) {
//   if (std::find(var_names.begin(), var_names.end(), var_name) ==
//       var_names.end()) {
//     var_names.push_back(var_name);
//   }
// }

// void ComplexEquationSystem::addTestVariableNameIfMissing(
//     std::string test_var_name) {
//   if (std::find(test_var_names.begin(), test_var_names.end(), test_var_name)
//   ==
//       test_var_names.end()) {
//     test_var_names.push_back(test_var_name);
//   }
// }

// void ComplexEquationSystem::addKernel(
//     std::string test_var_name,
//     hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel) {
//   addTestVariableNameIfMissing(test_var_name);
//   if (!blf_kernels_map.Has(test_var_name)) {
//     blf_kernels_map.Register(
//         test_var_name,
//         new mfem::Array<hephaestus::Kernel<mfem::ParBilinearForm> *>, true);
//   }
//   blf_kernels_map.Get(test_var_name)->Append(blf_kernel);
// }

// void ComplexEquationSystem::addKernel(
//     std::string test_var_name,
//     hephaestus::Kernel<mfem::ParSesquilinearForm> *lf_kernel) {
//   addTestVariableNameIfMissing(test_var_name);
//   if (!sqlf_kernels_map.Has(test_var_name)) {
//     sqlf_kernels_map.Register(
//         test_var_name,
//         new mfem::Array<hephaestus::Kernel<mfem::ParSesquilinearForm> *>,
//         true);
//   }

//   sqlf_kernels_map.Get(test_var_name)->Append(lf_kernel);
// }

// void ComplexEquationSystem::applyBoundaryConditions(hephaestus::BCMap
// &bc_map) {
//   ess_tdof_lists.resize(test_var_names.size());
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto test_var_name = test_var_names.at(i);
//     // Set default value of gridfunction used in essential BC. Values
//     // overwritten in applyEssentialBCs
//     *(xs.at(i)) = 0.0;
//     bc_map.applyEssentialBCs(test_var_name, ess_tdof_lists.at(i),
//     *(xs.at(i)),
//                              test_pfespaces.at(i)->GetParMesh());
//     bc_map.applyIntegratedBCs(test_var_name, *(sqlfs.Get(test_var_name)),
//                               test_pfespaces.at(i)->GetParMesh());
//   }
// }
// void ComplexEquationSystem::FormLinearSystem(mfem::OperatorHandle &op,
//                                              mfem::BlockVector &trueX,
//                                              mfem::BlockVector &trueRHS) {

//   // Allocate block operator
//   hBlocks.DeleteAll();
//   hBlocks.SetSize(test_var_names.size(), test_var_names.size());
//   // Form diagonal blocks.
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto &test_var_name = test_var_names.at(i);
//     auto blf = blfs.Get(test_var_name);
//     auto lf = sqlfs.Get(test_var_name);
//     mfem::Vector auxX, auxRHS;
//     hBlocks(i, i) = new mfem::HypreParMatrix;
//     blf->FormLinearSystem(ess_tdof_lists.at(i), *(xs.at(i)), *lf,
//                           *hBlocks(i, i), auxX, auxRHS);
//     trueX.GetBlock(i) = auxX;
//     trueRHS.GetBlock(i) = auxRHS;
//   }

//   // Form off-diagonal blocks
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto test_var_name = test_var_names.at(i);
//     for (int j = 0; j < test_var_names.size(); j++) {
//       auto trial_var_name = test_var_names.at(j);

//       mfem::Vector auxX, auxRHS;
//       mfem::ParSesquilinearForm auxLF(test_pfespaces.at(i));
//       auxLF = 0.0;
//       if (mblfs.Has(test_var_name) &&
//           mblfs.Get(test_var_name)->Has(trial_var_name)) {
//         auto mblf = mblfs.Get(test_var_name)->Get(trial_var_name);
//         hBlocks(i, j) = new mfem::HypreParMatrix;
//         mblf->FormRectangularLinearSystem(ess_tdof_lists.at(j),
//                                           ess_tdof_lists.at(i), *(xs.at(j)),
//                                           auxLF, *hBlocks(i, j), auxX,
//                                           auxRHS);
//         trueRHS.GetBlock(i) += auxRHS;
//       }
//     }
//   }
//   // Sync memory
//   for (int i = 0; i < test_var_names.size(); i++) {
//     trueX.GetBlock(0).SyncAliasMemory(trueX);
//     trueRHS.GetBlock(0).SyncAliasMemory(trueRHS);
//   }

//   // Create monolithic matrix
//   op.Reset(mfem::HypreParMatrixFromBlocks(hBlocks));
// }

// void ComplexEquationSystem::RecoverFEMSolution(
//     mfem::BlockVector &trueX,
//     mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto &test_var_name = test_var_names.at(i);
//     trueX.GetBlock(i).SyncAliasMemory(trueX);
//     variables.Get(test_var_name)->Distribute(&(trueX.GetBlock(i)));
//   }
// }

// void ComplexEquationSystem::Init(
//     mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//     const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//     hephaestus::BCMap &bc_map,
//     hephaestus::DomainProperties &domain_properties) {

//   // Add optional kernels to the ComplexEquationSystem
//   addKernels();

//   for (auto &test_var_name : test_var_names) {
//     if (!variables.Has(test_var_name)) {
//       MFEM_ABORT("Test variable "
//                  << test_var_name
//                  << " requested by equation system during initialisation was
//                  "
//                     "not found in variables");
//     }
//     // Store pointers to variable FESpaces
//     test_pfespaces.push_back(variables.Get(test_var_name)->ParFESpace());
//     // Create auxiliary gridfunctions for applying Dirichlet conditions
//     xs.push_back(
//         new
//         mfem::ParGridFunction(variables.Get(test_var_name)->ParFESpace()));
//   }

//   // Initialise bilinear forms

//   for (const auto &[test_var_name, blf_kernels] : blf_kernels_map.GetMap()) {
//     for (int i = 0; i < blf_kernels->Size(); i++) {
//       (*blf_kernels)[i]->Init(variables, fespaces, bc_map,
//       domain_properties);
//     }
//     blf_kernels->MakeDataOwner();
//   }
//   // Initialise linear form kernels
//   for (const auto &[test_var_name, sqlf_kernels] : sqlf_kernels_map.GetMap())
//   {
//     for (int i = 0; i < sqlf_kernels->Size(); i++) {
//       (*sqlf_kernels)[i]->Init(variables, fespaces, bc_map,
//       domain_properties);
//     }
//     sqlf_kernels->MakeDataOwner();
//   }
//   // Initialise nonlinear form kernels
//   for (const auto &[test_var_name, nlf_kernels] : nlf_kernels_map.GetMap()) {
//     for (int i = 0; i < nlf_kernels->Size(); i++) {
//       (*nlf_kernels)[i]->Init(variables, fespaces, bc_map,
//       domain_properties);
//     }
//     nlf_kernels->MakeDataOwner();
//   }
//   // Initialise mixed bilinear form kernels
//   for (const auto &[test_var_name, mblf_kernels_map] :
//        mblf_kernels_map_map.GetMap()) {
//     for (const auto &[trial_var_name, mblf_kernels] :
//          mblf_kernels_map->GetMap()) {
//       for (int i = 0; i < mblf_kernels->Size(); i++) {
//         (*mblf_kernels)[i]->Init(variables, fespaces, bc_map,
//                                  domain_properties);
//       }
//       mblf_kernels->MakeDataOwner();
//     }
//   }
// }

// void ComplexEquationSystem::buildLinearForms(hephaestus::BCMap &bc_map,
//                                              hephaestus::Sources &sources) {
//   // Register linear forms
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto test_var_name = test_var_names.at(i);
//     if (sqlfs.Has(test_var_name)) {
//       sqlfs.Deregister(test_var_name, true);
//     }
//     sqlfs.Register(test_var_name, new
//     mfem::ParSesquilinearForm(test_pfespaces.at(i)),
//                    true);
//     *(sqlfs.Get(test_var_name)) = 0.0;
//   }
//   // Apply boundary conditions
//   applyBoundaryConditions(bc_map);

//   for (auto &test_var_name : test_var_names) {
//     // Apply kernels
//     auto lf = sqlfs.Get(test_var_name);
//     auto sqlf_kernels = sqlf_kernels_map.Get(test_var_name);
//     if (sqlf_kernels != NULL) {
//       for (auto &lf_kernel : *sqlf_kernels) {
//         lf_kernel->Apply(lf);
//       }
//     }
//     // Assemble
//     lf->Assemble();
//     if (test_var_name == test_var_names.at(0)) {
//       sources.Apply(lf);
//     }
//   }
// }

// void ComplexEquationSystem::buildBilinearForms() {
//   // Register bilinear forms
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto test_var_name = test_var_names.at(i);
//     if (blfs.Has(test_var_name)) {
//       blfs.Deregister(test_var_name, true);
//     }
//     blfs.Register(test_var_name,
//                   new mfem::ParBilinearForm(test_pfespaces.at(i)), true);

//     // Apply kernels
//     auto blf = blfs.Get(test_var_name);
//     auto blf_kernels = blf_kernels_map.Get(test_var_name);
//     if (blf_kernels != NULL) {
//       for (auto &blf_kernel : *blf_kernels) {
//         blf_kernel->Apply(blf);
//       }
//     }
//     // Assemble
//     blf->Assemble();
//   }
// }

// void ComplexEquationSystem::buildMixedBilinearForms() {
//   // Register mixed linear forms. Note that not all combinations may
//   // have a kernel

//   // Create mblf for each test/trial pair
//   for (int i = 0; i < test_var_names.size(); i++) {
//     auto test_var_name = test_var_names.at(i);
//     if (mblfs.Has(test_var_name)) {
//       mblfs.Deregister(test_var_name, true);
//     }
//     mfem::NamedFieldsMap<mfem::ParMixedBilinearForm> *test_mblfs =
//         new mfem::NamedFieldsMap<mfem::ParMixedBilinearForm>;
//     for (int j = 0; j < test_var_names.size(); j++) {
//       auto trial_var_name = test_var_names.at(j);

//       // Register MixedBilinearForm if kernels exist for it, and assemble
//       // kernels
//       if (mblf_kernels_map_map.Has(test_var_name) &&
//           mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name)) {
//         auto mblf_kernels =
//             mblf_kernels_map_map.Get(test_var_name)->Get(trial_var_name);
//         mfem::ParMixedBilinearForm *mblf = new mfem::ParMixedBilinearForm(
//             test_pfespaces.at(j), test_pfespaces.at(i));
//         // Apply all mixed kernels with this test/trial pair
//         for (auto &mblf_kernel : *mblf_kernels) {
//           mblf_kernel->Apply(mblf);
//         }
//         // Assemble mixed bilinear forms
//         mblf->Assemble();
//         // Register mixed bilinear forms associated with a single trial
//         variable
//         // for the current test variable
//         test_mblfs->Register(trial_var_name, mblf, true);
//       }
//     }
//     // Register all mixed bilinear form sets associated with a single test
//     // variable
//     mblfs.Register(test_var_name, test_mblfs, true);
//   }
// }

// void ComplexEquationSystem::buildEquationSystem(hephaestus::BCMap &bc_map,
//                                                 hephaestus::Sources &sources)
//                                                 {
//   buildLinearForms(bc_map, sources);
//   buildBilinearForms();
//   buildMixedBilinearForms();
// }

// } // namespace hephaestus

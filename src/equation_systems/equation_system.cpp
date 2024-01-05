#include "equation_system.hpp"

namespace hephaestus {

EquationSystem::EquationSystem(const hephaestus::InputParameters &params)
    : var_names(), test_var_names(), test_pfespaces(), blfs(), lfs(), nlfs(),
      mblfs(), ess_tdof_lists(), xs() {}

EquationSystem::~EquationSystem() {
  hBlocks.DeleteAll();
  // blfs.DeleteData(true);
  // lfs.DeleteData(true);

  // for (const auto &[test_var_name, blf_kernels] : blf_kernels_map.GetMap()) {
  //   blf_kernels->DeleteAll();
  // }
  // blf_kernels_map.DeleteData(true);

  // for (const auto &[test_var_name, lf_kernels] : lf_kernels_map.GetMap()) {
  //   lf_kernels->DeleteAll();
  // }
  // lf_kernels_map.DeleteData(true);
}

void EquationSystem::addVariableNameIfMissing(std::string var_name) {
  if (std::find(var_names.begin(), var_names.end(), var_name) ==
      var_names.end()) {
    var_names.push_back(var_name);
  }
}

void EquationSystem::addTestVariableNameIfMissing(std::string test_var_name) {
  if (std::find(test_var_names.begin(), test_var_names.end(), test_var_name) ==
      test_var_names.end()) {
    test_var_names.push_back(test_var_name);
  }
}

void EquationSystem::addKernel(
    std::string test_var_name,
    hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel) {
  addTestVariableNameIfMissing(test_var_name);
  if (!blf_kernels_map.Has(test_var_name)) {
    blf_kernels_map.Register(
        test_var_name,
        new mfem::Array<hephaestus::Kernel<mfem::ParBilinearForm> *>, true);
  }
  blf_kernels_map.Get(test_var_name)->Append(blf_kernel);
}

void EquationSystem::addKernel(
    std::string test_var_name,
    hephaestus::Kernel<mfem::ParLinearForm> *lf_kernel) {
  addTestVariableNameIfMissing(test_var_name);
  if (!lf_kernels_map.Has(test_var_name)) {
    lf_kernels_map.Register(
        test_var_name,
        new mfem::Array<hephaestus::Kernel<mfem::ParLinearForm> *>, true);
  }

  lf_kernels_map.Get(test_var_name)->Append(lf_kernel);
}

void EquationSystem::addKernel(
    std::string test_var_name,
    hephaestus::Kernel<mfem::ParNonlinearForm> *nlf_kernel) {
  addTestVariableNameIfMissing(test_var_name);
  if (!nlf_kernels_map.Has(test_var_name)) {
    nlf_kernels_map.Register(
        test_var_name,
        new mfem::Array<hephaestus::Kernel<mfem::ParNonlinearForm> *>, true);
  }

  nlf_kernels_map.Get(test_var_name)->Append(nlf_kernel);
}

void EquationSystem::addKernel(
    std::string trial_var_name, std::string test_var_name,
    hephaestus::Kernel<mfem::ParMixedBilinearForm> *mblf_kernel) {
  addTestVariableNameIfMissing(test_var_name);
  // Register new mblf kernels map if not present for this test variable
  if (!mblf_kernels_map_map.Has(test_var_name)) {
    mblf_kernels_map_map.Register(
        test_var_name,
        new hephaestus::NamedFieldsMap<
            mfem::Array<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>>,
        true);
  }
  // Register new mblf kernels map if not present for the test/trial variable
  // pair
  if (!mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name)) {
    mblf_kernels_map_map.Get(test_var_name)
        ->Register(
            trial_var_name,
            new mfem::Array<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>,
            true);
  }

  mblf_kernels_map_map.Get(test_var_name)
      ->Get(trial_var_name)
      ->Append(mblf_kernel);
}

void EquationSystem::applyBoundaryConditions(hephaestus::BCMap &bc_map) {
  ess_tdof_lists.resize(test_var_names.size());
  for (int i = 0; i < test_var_names.size(); i++) {
    auto test_var_name = test_var_names.at(i);
    // Set default value of gridfunction used in essential BC. Values
    // overwritten in applyEssentialBCs
    *(xs.at(i)) = 0.0;
    bc_map.applyEssentialBCs(test_var_name, ess_tdof_lists.at(i), *(xs.at(i)),
                             test_pfespaces.at(i)->GetParMesh());
    bc_map.applyIntegratedBCs(test_var_name, *(lfs.Get(test_var_name)),
                              test_pfespaces.at(i)->GetParMesh());
  }
}
void EquationSystem::FormLinearSystem(mfem::OperatorHandle &op,
                                      mfem::BlockVector &trueX,
                                      mfem::BlockVector &trueRHS) {

  // Allocate block operator
  hBlocks.DeleteAll();
  hBlocks.SetSize(test_var_names.size(), test_var_names.size());
  // Form diagonal blocks.
  for (int i = 0; i < test_var_names.size(); i++) {
    auto &test_var_name = test_var_names.at(i);
    auto blf = blfs.Get(test_var_name);
    auto lf = lfs.Get(test_var_name);
    mfem::Vector auxX, auxRHS;
    hBlocks(i, i) = new mfem::HypreParMatrix;
    blf->FormLinearSystem(ess_tdof_lists.at(i), *(xs.at(i)), *lf,
                          *hBlocks(i, i), auxX, auxRHS);
    trueX.GetBlock(i) = auxX;
    trueRHS.GetBlock(i) = auxRHS;
  }

  // Form off-diagonal blocks
  for (int i = 0; i < test_var_names.size(); i++) {
    auto test_var_name = test_var_names.at(i);
    for (int j = 0; j < test_var_names.size(); j++) {
      auto trial_var_name = test_var_names.at(j);

      mfem::Vector auxX, auxRHS;
      mfem::ParLinearForm auxLF(test_pfespaces.at(i));
      auxLF = 0.0;
      if (mblfs.Has(test_var_name) &&
          mblfs.Get(test_var_name)->Has(trial_var_name)) {
        auto mblf = mblfs.Get(test_var_name)->Get(trial_var_name);
        hBlocks(i, j) = new mfem::HypreParMatrix;
        mblf->FormRectangularLinearSystem(ess_tdof_lists.at(j),
                                          ess_tdof_lists.at(i), *(xs.at(j)),
                                          auxLF, *hBlocks(i, j), auxX, auxRHS);
        trueRHS.GetBlock(i) += auxRHS;
      }
    }
  }
  // Sync memory
  for (int i = 0; i < test_var_names.size(); i++) {
    trueX.GetBlock(0).SyncAliasMemory(trueX);
    trueRHS.GetBlock(0).SyncAliasMemory(trueRHS);
  }

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(hBlocks));
}

void EquationSystem::RecoverFEMSolution(
    mfem::BlockVector &trueX, hephaestus::GridFunctions &gridfunctions) {
  for (int i = 0; i < test_var_names.size(); i++) {
    auto &test_var_name = test_var_names.at(i);
    trueX.GetBlock(i).SyncAliasMemory(trueX);
    gridfunctions.Get(test_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
}

void EquationSystem::Init(hephaestus::GridFunctions &gridfunctions,
                          const hephaestus::FESpaces &fespaces,
                          hephaestus::BCMap &bc_map,
                          hephaestus::Coefficients &coefficients) {

  // Add optional kernels to the EquationSystem
  addKernels();

  for (auto &test_var_name : test_var_names) {
    if (!gridfunctions.Has(test_var_name)) {
      MFEM_ABORT("Test variable "
                 << test_var_name
                 << " requested by equation system during initialisation was "
                    "not found in gridfunctions");
    }
    // Store pointers to variable FESpaces
    test_pfespaces.push_back(gridfunctions.Get(test_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for applying Dirichlet conditions
    xs.push_back(new mfem::ParGridFunction(
        gridfunctions.Get(test_var_name)->ParFESpace()));
  }

  // Initialise bilinear forms

  for (const auto &[test_var_name, blf_kernels] : blf_kernels_map.GetMap()) {
    for (int i = 0; i < blf_kernels->Size(); i++) {
      (*blf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
    }
    blf_kernels->MakeDataOwner();
  }
  // Initialise linear form kernels
  for (const auto &[test_var_name, lf_kernels] : lf_kernels_map.GetMap()) {
    for (int i = 0; i < lf_kernels->Size(); i++) {
      (*lf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
    }
    lf_kernels->MakeDataOwner();
  }
  // Initialise nonlinear form kernels
  for (const auto &[test_var_name, nlf_kernels] : nlf_kernels_map.GetMap()) {
    for (int i = 0; i < nlf_kernels->Size(); i++) {
      (*nlf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
    }
    nlf_kernels->MakeDataOwner();
  }
  // Initialise mixed bilinear form kernels
  for (const auto &[test_var_name, mblf_kernels_map] :
       mblf_kernels_map_map.GetMap()) {
    for (const auto &[trial_var_name, mblf_kernels] :
         mblf_kernels_map->GetMap()) {
      for (int i = 0; i < mblf_kernels->Size(); i++) {
        (*mblf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
      }
      mblf_kernels->MakeDataOwner();
    }
  }
}

void EquationSystem::buildLinearForms(hephaestus::BCMap &bc_map,
                                      hephaestus::Sources &sources) {
  // Register linear forms
  for (int i = 0; i < test_var_names.size(); i++) {
    auto test_var_name = test_var_names.at(i);
    if (lfs.Has(test_var_name)) {
      lfs.Deregister(test_var_name);
    }
    lfs.Register(test_var_name, new mfem::ParLinearForm(test_pfespaces.at(i)),
                 true);
    *(lfs.Get(test_var_name)) = 0.0;
  }
  // Apply boundary conditions
  applyBoundaryConditions(bc_map);

  for (auto &test_var_name : test_var_names) {
    // Apply kernels
    auto lf = lfs.Get(test_var_name);
    // Assemble. Must be done before applying kernels that add to lf.
    lf->Assemble();

    auto lf_kernels = lf_kernels_map.Get(test_var_name);
    if (lf_kernels != NULL) {
      for (auto &lf_kernel : *lf_kernels) {
        lf_kernel->Apply(lf);
      }
    }
    if (test_var_name == test_var_names.at(0)) {
      sources.Apply(lf);
    }
  }
}

void EquationSystem::buildBilinearForms() {
  // Register bilinear forms
  for (int i = 0; i < test_var_names.size(); i++) {
    auto test_var_name = test_var_names.at(i);
    if (blfs.Has(test_var_name)) {
      blfs.Deregister(test_var_name);
    }
    blfs.Register(test_var_name,
                  new mfem::ParBilinearForm(test_pfespaces.at(i)), true);

    // Apply kernels
    auto blf = blfs.Get(test_var_name);
    auto blf_kernels = blf_kernels_map.Get(test_var_name);
    if (blf_kernels != NULL) {
      for (auto &blf_kernel : *blf_kernels) {
        blf_kernel->Apply(blf);
      }
    }
    // Assemble
    blf->Assemble();
  }
}

void EquationSystem::buildMixedBilinearForms() {
  // Register mixed linear forms. Note that not all combinations may
  // have a kernel

  // Create mblf for each test/trial pair
  for (int i = 0; i < test_var_names.size(); i++) {
    auto test_var_name = test_var_names.at(i);
    if (mblfs.Has(test_var_name)) {
      mblfs.Deregister(test_var_name);
    }
    hephaestus::NamedFieldsMap<mfem::ParMixedBilinearForm> *test_mblfs =
        new hephaestus::NamedFieldsMap<mfem::ParMixedBilinearForm>;
    for (int j = 0; j < test_var_names.size(); j++) {
      auto trial_var_name = test_var_names.at(j);

      // Register MixedBilinearForm if kernels exist for it, and assemble
      // kernels
      if (mblf_kernels_map_map.Has(test_var_name) &&
          mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name)) {
        auto mblf_kernels =
            mblf_kernels_map_map.Get(test_var_name)->Get(trial_var_name);
        mfem::ParMixedBilinearForm *mblf = new mfem::ParMixedBilinearForm(
            test_pfespaces.at(j), test_pfespaces.at(i));
        // Apply all mixed kernels with this test/trial pair
        for (auto &mblf_kernel : *mblf_kernels) {
          mblf_kernel->Apply(mblf);
        }
        // Assemble mixed bilinear forms
        mblf->Assemble();
        // Register mixed bilinear forms associated with a single trial variable
        // for the current test variable
        test_mblfs->Register(trial_var_name, mblf, true);
      }
    }
    // Register all mixed bilinear form sets associated with a single test
    // variable
    mblfs.Register(test_var_name, test_mblfs, true);
  }
}

void EquationSystem::buildEquationSystem(hephaestus::BCMap &bc_map,
                                         hephaestus::Sources &sources) {
  buildLinearForms(bc_map, sources);
  buildBilinearForms();
  buildMixedBilinearForms();
}

TimeDependentEquationSystem::TimeDependentEquationSystem(
    const hephaestus::InputParameters &params)
    : EquationSystem(params), var_time_derivative_names(), dtCoef(1.0) {}

void TimeDependentEquationSystem::addVariableNameIfMissing(
    std::string var_name) {
  EquationSystem::addVariableNameIfMissing(var_name);
  std::string var_time_derivative_name = GetTimeDerivativeName(var_name);
  if (std::find(var_time_derivative_names.begin(),
                var_time_derivative_names.end(),
                var_time_derivative_name) == var_time_derivative_names.end()) {
    var_time_derivative_names.push_back(var_time_derivative_name);
  }
}

void TimeDependentEquationSystem::setTimeStep(double dt) {
  if (fabs(dt - dtCoef.constant) > 1.0e-12 * dt) {
    dtCoef.constant = dt;
    for (int i = 0; i < test_var_names.size(); i++) {
      auto test_var_name = test_var_names.at(i);
      auto blf = blfs.Get(test_var_name);
      blf->Update();
      blf->Assemble();
    }
  }
}

void TimeDependentEquationSystem::updateEquationSystem(
    hephaestus::BCMap &bc_map, hephaestus::Sources &sources) {
  buildLinearForms(bc_map, sources);
};

} // namespace hephaestus

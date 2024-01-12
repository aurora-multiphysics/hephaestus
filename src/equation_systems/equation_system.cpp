#include "equation_system.hpp"

namespace hephaestus {

EquationSystem::EquationSystem(const hephaestus::InputParameters &params)
    : var_names(), test_var_names(), test_pfespaces(), blfs(), lfs(), nlfs(),
      mblfs(), ess_tdof_lists(), xs() {}

EquationSystem::~EquationSystem() {
  fprintf(stdout, "Calling %s\n", __func__);
  fflush(stdout);
  hBlocks.DeleteAll();

  for (const auto par_grid_function : xs) {
    delete par_grid_function;
  }
}

bool EquationSystem::vectorContainsName(
    const std::vector<std::string> &the_vector, const std::string &name) const {

  auto iter = std::find(the_vector.begin(), the_vector.end(), name);

  return (iter != the_vector.end());
}

void EquationSystem::addVariableNameIfMissing(const std::string &var_name) {
  if (!vectorContainsName(var_names, var_name)) {
    var_names.push_back(var_name);
  }
}

void EquationSystem::addTestVariableNameIfMissing(
    const std::string &test_var_name) {
  if (!vectorContainsName(test_var_names, test_var_name)) {
    test_var_names.push_back(test_var_name);
  }
}

void EquationSystem::addKernel(
    const std::string &test_var_name,
    std::unique_ptr<ParBilinearFormKernel> blf_kernel) {

  addTestVariableNameIfMissing(test_var_name);

  if (!blf_kernels_map.Has(test_var_name)) {
    // 1. Create kernels vector.
    auto kernels = new std::vector<decltype(blf_kernel)>;

    // 2. Register with map to prevent leaks.
    blf_kernels_map.Register(test_var_name, kernels, true);
  }

  blf_kernels_map.Get(test_var_name)->push_back(std::move(blf_kernel));
}

void EquationSystem::addKernel(const std::string &test_var_name,
                               std::unique_ptr<ParLinearFormKernel> lf_kernel) {

  addTestVariableNameIfMissing(test_var_name);

  if (!lf_kernels_map.Has(test_var_name)) {
    auto kernels = new std::vector<decltype(lf_kernel)>;

    lf_kernels_map.Register(test_var_name, kernels, true);
  }

  lf_kernels_map.Get(test_var_name)->push_back(std::move(lf_kernel));
}

void EquationSystem::addKernel(
    const std::string &test_var_name,
    std::unique_ptr<ParNonlinearFormKernel> nlf_kernel) {

  addTestVariableNameIfMissing(test_var_name);

  if (!nlf_kernels_map.Has(test_var_name)) {
    auto kernels = new std::vector<decltype(nlf_kernel)>;

    nlf_kernels_map.Register(test_var_name, kernels, true);
  }

  nlf_kernels_map.Get(test_var_name)->push_back(std::move(nlf_kernel));
}

void EquationSystem::addKernel(
    const std::string &trial_var_name, const std::string &test_var_name,
    std::unique_ptr<ParMixedBilinearFormKernel> mblf_kernel) {

  addTestVariableNameIfMissing(test_var_name);

  // Register new mblf kernels map if not present for this test variable
  if (!mblf_kernels_map_map.Has(test_var_name)) {
    auto kernel_field_map =
        new hephaestus::NamedFieldsMap<std::vector<decltype(mblf_kernel)>>;

    mblf_kernels_map_map.Register(test_var_name, kernel_field_map, true);
  }

  // Register new mblf kernels map if not present for the test/trial variable
  // pair
  if (!mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name)) {
    auto kernels = new std::vector<std::unique_ptr<ParMixedBilinearFormKernel>>;

    mblf_kernels_map_map.Get(test_var_name)
        ->Register(trial_var_name, kernels, true);
  }

  mblf_kernels_map_map.Get(test_var_name)
      ->Get(trial_var_name)
      ->push_back(std::move(mblf_kernel));
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
    for (int i = 0; i < blf_kernels->size(); i++) {
      (*blf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
    }
  }
  // Initialise linear form kernels
  for (const auto &[test_var_name, lf_kernels] : lf_kernels_map.GetMap()) {
    for (int i = 0; i < lf_kernels->size(); i++) {
      (*lf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
    }
  }
  // Initialise nonlinear form kernels
  for (const auto &[test_var_name, nlf_kernels] : nlf_kernels_map.GetMap()) {
    for (int i = 0; i < nlf_kernels->size(); i++) {
      (*nlf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
    }
  }
  // Initialise mixed bilinear form kernels
  for (const auto &[test_var_name, mblf_kernels_map] :
       mblf_kernels_map_map.GetMap()) {
    for (const auto &[trial_var_name, mblf_kernels] :
         mblf_kernels_map->GetMap()) {
      for (int i = 0; i < mblf_kernels->size(); i++) {
        (*mblf_kernels)[i]->Init(gridfunctions, fespaces, bc_map, coefficients);
      }
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
    if (lf_kernels != nullptr) {
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
    if (blf_kernels != nullptr) {
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
    const std::string &var_name) {
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
  buildBilinearForms();
  buildMixedBilinearForms();
};

} // namespace hephaestus

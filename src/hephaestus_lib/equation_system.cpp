#include "equation_system.hpp"

namespace hephaestus {

EquationSystem::EquationSystem(const hephaestus::InputParameters &params)
    : var_names(params.GetParam<std::vector<std::string>>("VariableNames")),
      test_var_names(var_names), test_pfespaces(), blfs(), lfs(), nlfs(),
      mblfs(), ess_tdof_lists(test_var_names.size()), xs() {}

void EquationSystem::addKernel(
    std::string test_var_name,
    hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel) {
  if (!blf_kernels_map.Has(test_var_name)) {
    blf_kernels_map.Register(
        test_var_name,
        new std::vector<hephaestus::Kernel<mfem::ParBilinearForm> *>, true);
  }
  blf_kernels_map.Get(test_var_name)->push_back(blf_kernel);
}

void EquationSystem::addKernel(
    std::string test_var_name,
    hephaestus::Kernel<mfem::ParLinearForm> *lf_kernel) {
  if (!lf_kernels_map.Has(test_var_name)) {
    lf_kernels_map.Register(
        test_var_name,
        new std::vector<hephaestus::Kernel<mfem::ParLinearForm> *>, true);
  }

  lf_kernels_map.Get(test_var_name)->push_back(lf_kernel);
}

void EquationSystem::addKernel(
    std::string test_var_name,
    hephaestus::Kernel<mfem::ParNonlinearForm> *nlf_kernel) {
  if (!nlf_kernels_map.Has(test_var_name)) {
    nlf_kernels_map.Register(
        test_var_name,
        new std::vector<hephaestus::Kernel<mfem::ParNonlinearForm> *>, true);
  }

  nlf_kernels_map.Get(test_var_name)->push_back(nlf_kernel);
}

void EquationSystem::addKernel(
    std::string trial_var_name, std::string test_var_name,
    hephaestus::Kernel<mfem::ParMixedBilinearForm> *mblf_kernel) {
  // Register new mblf kernels map if not present for this test variable
  if (!mblf_kernels_map_map.Has(test_var_name)) {
    mblf_kernels_map_map.Register(
        test_var_name,
        new mfem::NamedFieldsMap<
            std::vector<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>>,
        true);
  }
  // Register new mblf kernels map if not present for the test/trial variable
  // pair
  if (!mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name)) {
    mblf_kernels_map_map.Get(test_var_name)
        ->Register(
            trial_var_name,
            new std::vector<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>,
            true);
  }

  mblf_kernels_map_map.Get(test_var_name)
      ->Get(trial_var_name)
      ->push_back(mblf_kernel);
}

void EquationSystem::applyBoundaryConditions(hephaestus::BCMap &bc_map) {
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
  hBlocks.SetSize(test_var_names.size(), test_var_names.size());
  hBlocks = NULL;
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

void EquationSystem::RecoverFEMSolution(mfem::Vector &X,
                                        mfem::ParGridFunction &test_variable) {
  for (int i = 0; i < test_var_names.size(); i++) {
    auto &test_var_name = test_var_names.at(i);
    auto blf = blfs.Get(test_var_name);
    auto lf = lfs.Get(test_var_name);
    blf->RecoverFEMSolution(X, *lf, test_variable);
  }
}

void EquationSystem::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  for (auto &test_var_name : test_var_names) {
    // Store pointers to variable FESpaces
    test_pfespaces.push_back(variables.Get(test_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for applying Dirichlet conditions
    xs.push_back(
        new mfem::ParGridFunction(variables.Get(test_var_name)->ParFESpace()));
  }
  // Add optional kernels to the EquationSystem
  addKernels();

  // for (auto &test_var_name : test_var_names) {

  // Initialise bilinear forms
  for (const auto &[test_var_name, blf_kernels] : blf_kernels_map.GetMap()) {
    for (int i = 0; i < blf_kernels->size(); i++) {
      blf_kernels->at(i)->Init(variables, fespaces, bc_map, domain_properties);
    }
  }
  // Initialise linear forms
  for (const auto &[test_var_name, lf_kernels] : lf_kernels_map.GetMap()) {
    for (int i = 0; i < lf_kernels->size(); i++) {
      lf_kernels->at(i)->Init(variables, fespaces, bc_map, domain_properties);
    }
  }
  // Initialise nonlinear forms
  for (const auto &[test_var_name, nlf_kernels] : nlf_kernels_map.GetMap()) {
    for (int i = 0; i < nlf_kernels->size(); i++) {
      nlf_kernels->at(i)->Init(variables, fespaces, bc_map, domain_properties);
    }
  }
  // Initialise mixed bilinear forms
  for (const auto &[test_var_name, mblf_kernels_map] :
       mblf_kernels_map_map.GetMap()) {
    for (const auto &[trial_var_name, mblf_kernels] :
         mblf_kernels_map->GetMap()) {
      for (int i = 0; i < mblf_kernels->size(); i++) {
        mblf_kernels->at(i)->Init(variables, fespaces, bc_map,
                                  domain_properties);
      }
    }
  }
  // }
}

void EquationSystem::buildLinearForms(hephaestus::BCMap &bc_map,
                                      hephaestus::Sources &sources) {
  // Register linear forms
  for (int i = 0; i < test_var_names.size(); i++) {
    auto test_var_name = test_var_names.at(i);
    if (lfs.Has(test_var_name)) {
      lfs.Deregister(test_var_name, true);
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
    auto lf_kernels = lf_kernels_map.Get(test_var_name);
    if (lf_kernels != NULL) {
      for (auto &lf_kernel : *lf_kernels) {
        lf_kernel->Apply(lf);
      }
    }
    // Assemble
    lf->Assemble();
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
      blfs.Deregister(test_var_name, true);
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
      mblfs.Deregister(test_var_name, true);
    }
    mfem::NamedFieldsMap<mfem::ParMixedBilinearForm> *test_mblfs =
        new mfem::NamedFieldsMap<mfem::ParMixedBilinearForm>;
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
    : EquationSystem(params),
      var_time_derivative_names(
          params.GetParam<std::vector<std::string>>("TimeDerivativeNames")),
      dtCoef(1.0) {}

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

CurlCurlEquationSystem::CurlCurlEquationSystem(
    const hephaestus::InputParameters &params)
    : TimeDependentEquationSystem(params),
      alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
      beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
      dtalpha_coef_name(std::string("dt_") + alpha_coef_name) {

  test_var_names.at(0) = var_time_derivative_names.at(0);
}

void CurlCurlEquationSystem::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  domain_properties.scalar_property_map[dtalpha_coef_name] =
      new mfem::TransformedCoefficient(
          &dtCoef, domain_properties.scalar_property_map[alpha_coef_name],
          prodFunc);
  TimeDependentEquationSystem::Init(variables, fespaces, bc_map,
                                    domain_properties);
}

void CurlCurlEquationSystem::addKernels() {
  // (α∇×u_{n}, ∇×u')
  hephaestus::InputParameters weakCurlCurlParams;
  weakCurlCurlParams.SetParam("VariableName", var_names.at(0));
  weakCurlCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(test_var_names.at(0),
            new hephaestus::WeakCurlCurlKernel(weakCurlCurlParams));

  // (αdt∇×du/dt_{n+1}, ∇×u')
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("VariableName", test_var_names.at(0));
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(test_var_names.at(0),
            new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βdu/dt_{n+1}, u')
  hephaestus::InputParameters vectorFEMassParams;
  vectorFEMassParams.SetParam("VariableName", test_var_names.at(0));
  vectorFEMassParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(test_var_names.at(0),
            new hephaestus::VectorFEMassKernel(vectorFEMassParams));
}

AVEquationSystem::AVEquationSystem(const hephaestus::InputParameters &params)
    : TimeDependentEquationSystem(params),
      alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
      beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
      dtalpha_coef_name(std::string("dt_") + alpha_coef_name),
      neg_beta_coef_name(std::string("negative_") + beta_coef_name),
      negCoef(-1.0) {
  test_var_names.at(0) = var_time_derivative_names.at(0);
  test_var_names.at(1) = var_names.at(1);
}

void AVEquationSystem::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  domain_properties.scalar_property_map[dtalpha_coef_name] =
      new mfem::TransformedCoefficient(
          &dtCoef, domain_properties.scalar_property_map[alpha_coef_name],
          prodFunc);
  domain_properties.scalar_property_map[neg_beta_coef_name] =
      new mfem::TransformedCoefficient(
          &negCoef, domain_properties.scalar_property_map[beta_coef_name],
          prodFunc);
  TimeDependentEquationSystem::Init(variables, fespaces, bc_map,
                                    domain_properties);
}

void AVEquationSystem::addKernels() {
  // (α∇×A_{n}, ∇×A') - careful about trial var name!
  hephaestus::InputParameters weakCurlCurlParams;
  weakCurlCurlParams.SetParam("VariableName", var_names.at(0));
  weakCurlCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(test_var_names.at(0),
            new hephaestus::WeakCurlCurlKernel(weakCurlCurlParams));

  // (αdt∇×dA/dt_{n+1}, ∇×A')
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("VariableName", test_var_names.at(0));
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(test_var_names.at(0),
            new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βdA/dt_{n+1}, A')
  hephaestus::InputParameters vectorFEMassParams;
  vectorFEMassParams.SetParam("VariableName", test_var_names.at(0));
  vectorFEMassParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(test_var_names.at(0),
            new hephaestus::VectorFEMassKernel(vectorFEMassParams));

  // (σ ∇ V, dA'/dt)
  hephaestus::InputParameters mixedVectorGradientParams;
  mixedVectorGradientParams.SetParam("VariableName", test_var_names.at(0));
  mixedVectorGradientParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(
      test_var_names.at(1), test_var_names.at(0),
      new hephaestus::MixedVectorGradientKernel(mixedVectorGradientParams));

  // (σ ∇ V, ∇ V')
  hephaestus::InputParameters diffusionParams;
  diffusionParams.SetParam("VariableName", test_var_names.at(1));
  diffusionParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(test_var_names.at(1),
            new hephaestus::DiffusionKernel(diffusionParams));

  // (σdA/dt, ∇ V')
  hephaestus::InputParameters vectorFEWeakDivergenceParams;
  vectorFEWeakDivergenceParams.SetParam("VariableName", test_var_names.at(1));
  vectorFEWeakDivergenceParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(test_var_names.at(0), test_var_names.at(1),
            new hephaestus::VectorFEWeakDivergenceKernel(
                vectorFEWeakDivergenceParams));
}

} // namespace hephaestus

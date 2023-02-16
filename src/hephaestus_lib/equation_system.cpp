#include "equation_system.hpp"

namespace hephaestus {

WeakForm::WeakForm(const hephaestus::InputParameters &params)
    : _test_var_name(params.GetParam<std::string>("TestVariableName")),
      test_pfes(NULL), blf(NULL), lf(NULL), nlf(NULL), mblfs(), ess_tdof_list(),
      x(NULL) {}

void WeakForm::addKernel(
    hephaestus::Kernel<mfem::ParBilinearForm> *blf_kernel) {
  blf_kernels.Append(blf_kernel);
}

void WeakForm::addKernel(hephaestus::Kernel<mfem::ParLinearForm> *lf_kernel) {
  lf_kernels.Append(lf_kernel);
}

void WeakForm::addKernel(
    hephaestus::Kernel<mfem::ParNonlinearForm> *nlf_kernel) {
  nlf_kernels.Append(nlf_kernel);
}

void WeakForm::addKernel(
    std::string trial_var_name,
    hephaestus::Kernel<mfem::ParMixedBilinearForm> *mblf_kernel) {
  if (!mblf_kernels_map.Has(trial_var_name)) {
    mblf_kernels_map.Register(
        trial_var_name,
        new mfem::Array<hephaestus::Kernel<mfem::ParMixedBilinearForm> *>,
        true);
  }
  mblf_kernels_map.Get(trial_var_name)->Append(mblf_kernel);
}

void WeakForm::applyBoundaryConditions(hephaestus::BCMap &bc_map) {
  *x = 0.0;
  bc_map.applyEssentialBCs(_test_var_name, ess_tdof_list, *x,
                           test_pfes->GetParMesh());
  bc_map.applyIntegratedBCs(_test_var_name, *lf, test_pfes->GetParMesh());
}

void WeakForm::FormLinearSystem(mfem::HypreParMatrix &A, mfem::Vector &X,
                                mfem::Vector &B) {
  blf->FormLinearSystem(ess_tdof_list, *x, *lf, A, X, B);
}

void WeakForm::RecoverFEMSolution(mfem::Vector &X,
                                  mfem::ParGridFunction &test_variable) {
  blf->RecoverFEMSolution(X, *lf, test_variable);
}

void WeakForm::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  test_pfes = variables.Get(_test_var_name)->ParFESpace();
  x = new mfem::ParGridFunction(test_pfes);
  addKernels(variables, fespaces, bc_map, domain_properties);

  for (int i = 0; i < blf_kernels.Size(); i++) {
    blf_kernels[i]->Init(variables, fespaces, bc_map, domain_properties);
  }
  for (int i = 0; i < lf_kernels.Size(); i++) {
    lf_kernels[i]->Init(variables, fespaces, bc_map, domain_properties);
  }
  for (int i = 0; i < nlf_kernels.Size(); i++) {
    nlf_kernels[i]->Init(variables, fespaces, bc_map, domain_properties);
  }
  for (const auto &[trial_var_name, mblf_kernels] : mblf_kernels_map.GetMap()) {
    for (int i = 0; i < mblf_kernels->Size(); i++) {
      (*mblf_kernels)[i]->Init(variables, fespaces, bc_map, domain_properties);
    }
  }
}

void WeakForm::buildLinearForm(hephaestus::BCMap &bc_map,
                               hephaestus::Sources &sources) {
  if (lf != NULL) {
    delete lf;
  }
  lf = new mfem::ParLinearForm(test_pfes);
  applyBoundaryConditions(bc_map);
  for (int i = 0; i < lf_kernels.Size(); i++) {
    lf_kernels[i]->Apply(lf);
  }
  lf->Assemble();
  sources.Apply(lf);
}

void WeakForm::buildBilinearForm() {
  if (blf != NULL) {
    delete blf;
  }
  blf = new mfem::ParBilinearForm(test_pfes);
  for (int i = 0; i < blf_kernels.Size(); i++) {
    blf_kernels[i]->Apply(blf);
  }
  blf->Assemble();
}

void WeakForm::buildWeakForm(hephaestus::BCMap &bc_map,
                             hephaestus::Sources &sources) {
  buildLinearForm(bc_map, sources);
  buildBilinearForm();
}

TimeDependentWeakForm::TimeDependentWeakForm(
    const hephaestus::InputParameters &params)
    : WeakForm(params), dtCoef(1.0) {}

void TimeDependentWeakForm::setTimeStep(double dt) { dtCoef.constant = dt; }

void TimeDependentWeakForm::updateWeakForm(hephaestus::BCMap &bc_map,
                                           hephaestus::Sources &sources) {
  buildLinearForm(bc_map, sources);
};

CurlCurlWeakForm::CurlCurlWeakForm(const hephaestus::InputParameters &params)
    : TimeDependentWeakForm(params),
      alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
      beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
      dtalpha_coef_name(std::string("dt_") + alpha_coef_name) {}

void CurlCurlWeakForm::addKernels(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {

  domain_properties.scalar_property_map[dtalpha_coef_name] =
      new mfem::TransformedCoefficient(
          &dtCoef, domain_properties.scalar_property_map[alpha_coef_name],
          prodFunc);

  // (α∇×u_{n}, ∇×u')
  hephaestus::InputParameters weakCurlCurlParams;
  weakCurlCurlParams.SetParam("VariableName", _test_var_name);
  weakCurlCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(new hephaestus::WeakCurlCurlKernel(weakCurlCurlParams));

  // (αdt∇×du/dt_{n+1}, ∇×u')
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("VariableName", _test_var_name);
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βdu/dt_{n+1}, u')
  hephaestus::InputParameters vectorFEMassParams;
  vectorFEMassParams.SetParam("VariableName", _test_var_name);
  vectorFEMassParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(new hephaestus::VectorFEMassKernel(vectorFEMassParams));
}

void CurlCurlWeakForm::updateWeakForm(hephaestus::BCMap &bc_map,
                                      hephaestus::Sources &sources) {
  buildLinearForm(bc_map, sources);
};

void CurlCurlWeakForm::setTimeStep(double dt) {
  if (blf == NULL || fabs(dt - dtCoef.constant) > 1.0e-12 * dt) {
    TimeDependentWeakForm::setTimeStep(dt);
    blf->Update();
    blf->Assemble();
  }
}
// EquationSystem::EquationSystem(NamedFieldsMap<hephaestus::WeakForm> eqns) {

// }

// EquationSystem::EquationSystem(eqn_names, variables) {
//   for (const auto &[var_name, variables] : GetMap()) {
//   }

//   for (eqn_name in eqn_names) {
//     eqn = factory_buildEquation(eqn_name.Get(test_var_name));
//     equations.Register(
//         eqn_name, test_var_name,
//         new mfem::ParMixedBilinearForm(trial_variable.ParFESpace(),
//         test_pfes), false);

//     factory.buildEquation(eqn_name)
//   }

//   true_offsets.SetSize(vars.Size());
//   true_offsets[0] = 0;
//   i = 1;
//   for (test_var in vars) {
//     true_offsets[i] = test_var->ParFESpace()->GetVSize();
//     i += 1;
//   }
//   true_offsets.PartialSum();

//   this->height = true_offsets[vars.Size() - 1];
//   this->width = true_offsets[vars.Size() - 1];
// }

// // CurlCurlWeakForm::buildMixedBilinearForm(const std::string
// trial_var_name,
// //                                       mfem::ParGridFunction
// &trial_variable)
// //                                       {
// //   mblfs.Register(
// //       trial_var_name,
// //       new mfem::ParMixedBilinearForm(trial_variable.ParFESpace(),
// test_pfes),
// //       false);
// //   mblf = mblfs.Get(trial_var_name);
// //   // act on mblf
// // }

} // namespace hephaestus

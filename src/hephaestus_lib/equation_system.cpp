#include "equation_system.hpp"

namespace hephaestus {

WeakForm::WeakForm(const std::string test_var_name,
                   mfem::ParGridFunction &test_variable)
    : _test_var_name(test_var_name), test_pfes(test_variable.ParFESpace()),
      blf(new mfem::ParBilinearForm(test_pfes)),
      lf(new mfem::ParLinearForm(test_pfes)), nlf(NULL), mblfs(),
      ess_tdof_list(), x(test_pfes) {
  *lf = 0.0;
}

void WeakForm::applyBoundaryConditions(hephaestus::BCMap &bc_map) {
  x = 0.0;
  bc_map.applyEssentialBCs(_test_var_name, ess_tdof_list, x,
                           test_pfes->GetParMesh());
  bc_map.applyIntegratedBCs(_test_var_name, *lf, test_pfes->GetParMesh());
}

void WeakForm::FormLinearSystem(mfem::HypreParMatrix &A, mfem::Vector &X,
                                mfem::Vector &B) {
  blf->FormLinearSystem(ess_tdof_list, x, *lf, A, X, B);
}

void WeakForm::RecoverFEMSolution(mfem::Vector &X,
                                  mfem::ParGridFunction &test_variable) {
  blf->RecoverFEMSolution(X, *lf, test_variable);
}

TimeDependentWeakForm::TimeDependentWeakForm(
    const std::string test_var_name, mfem::ParGridFunction &test_variable)
    : WeakForm(test_var_name, test_variable), dtCoef(1.0) {}

void TimeDependentWeakForm::setTimeStep(double dt) { dtCoef.constant = dt; }

void TimeDependentWeakForm::updateWeakForm(hephaestus::BCMap &bc_map,
                                           hephaestus::Sources &sources) {
  buildLinearForm(bc_map, sources);
};

CurlCurlWeakForm::CurlCurlWeakForm(const std::string test_var_name,
                                   mfem::ParGridFunction &test_variable,
                                   mfem::ParGridFunction &coupled_variable,
                                   mfem::Coefficient *alphaCoef_,
                                   mfem::Coefficient *betaCoef_)
    : TimeDependentWeakForm(test_var_name, test_variable),
      curlCurl(new mfem::ParBilinearForm(test_pfes)), u_(coupled_variable),
      alphaCoef(alphaCoef_), betaCoef(betaCoef_),
      dtAlphaCoef(
          new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc)) {}

void CurlCurlWeakForm::buildLinearForm(hephaestus::BCMap &bc_map,
                                       hephaestus::Sources &sources) {
  if (curlCurl != NULL) {
    delete curlCurl;
  }
  curlCurl = new mfem::ParBilinearForm(test_pfes);
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*alphaCoef));
  curlCurl->Assemble();

  if (lf != NULL) {
    delete lf;
  }
  lf = new mfem::ParLinearForm(test_pfes);
  applyBoundaryConditions(bc_map);
  lf->Assemble();

  curlCurl->AddMultTranspose(u_, *lf, -1.0);
  sources.Apply(lf);
}

void CurlCurlWeakForm::buildBilinearForm() {
  if (blf != NULL) {
    delete blf;
  }
  blf = new mfem::ParBilinearForm(test_pfes);
  blf->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*betaCoef));
  blf->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*dtAlphaCoef));
  blf->Assemble();
}

void CurlCurlWeakForm::buildWeakForm(hephaestus::BCMap &bc_map,
                                     hephaestus::Sources &sources) {
  buildLinearForm(bc_map, sources);
  buildBilinearForm();
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

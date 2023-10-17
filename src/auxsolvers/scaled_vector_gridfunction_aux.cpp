#include "scaled_vector_gridfunction_aux.hpp"

namespace hephaestus {

ScaledVectorGridFunctionAux::ScaledVectorGridFunctionAux(
    const hephaestus::InputParameters &params)
    : AuxSolver(), coef(nullptr), input_gf(nullptr), scaled_gf(nullptr),
      coef_name(params.GetParam<std::string>("CoefficientName")),
      input_gf_name(params.GetParam<std::string>("InputVariableName")),
      scaled_gf_name(params.GetParam<std::string>("ScaledVariableName")),
      solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters())),
      test_fes(nullptr), trial_fes(nullptr), a(nullptr), a_mixed(nullptr),
      a_mat(nullptr), mixed_mat(nullptr), solver(nullptr) {}

void ScaledVectorGridFunctionAux::Init(
    const hephaestus::GridFunctions &gridfunctions,
    hephaestus::Coefficients &coefficients) {
  input_gf = gridfunctions.Get(input_gf_name);
  if (input_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << input_gf_name
               << " not found when initializing ScaledVectorGridFunctionAux");
  }
  scaled_gf = gridfunctions.Get(scaled_gf_name);
  if (scaled_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << scaled_gf_name
               << " not found when initializing ScaledVectorGridFunctionAux");
  }
  coef = dynamic_cast<mfem::Coefficient *>(coefficients.scalars.Get(coef_name));
  if (coef == NULL) {
    MFEM_ABORT("Coefficient "
               << coef_name
               << " not found when initializing ScaledVectorGridFunctionAux");
  }

  test_fes = scaled_gf->ParFESpace();
  trial_fes = input_gf->ParFESpace();
  mixed_mat = a_mixed->ParallelAssemble();
  a_mat = a->ParallelAssemble();
  solver = new hephaestus::DefaultH1PCGSolver(solver_options, *a_mat);
}

void ScaledVectorGridFunctionAux::buildBilinearForm() {
  a = new mfem::ParBilinearForm(test_fes);
  a->AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  a->Assemble();
}

void ScaledVectorGridFunctionAux::buildMixedBilinearForm() {
  a_mixed = new mfem::ParMixedBilinearForm(trial_fes, test_fes);
  a_mixed->AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(*coef));
  a_mixed->Assemble();
}

void ScaledVectorGridFunctionAux::Solve(double t) {
  mfem::Vector B(test_fes->GetTrueVSize());  // Linear form true DOFs
  mfem::Vector X(test_fes->GetTrueVSize());  // H(Div) gridfunction true DOFs
  mfem::Vector P(trial_fes->GetTrueVSize()); // H(Curl) gridfunction true DOFs
  input_gf->GetTrueDofs(P);
  mixed_mat->Mult(P, B);
  solver->Mult(B, X);
  scaled_gf->SetFromTrueDofs(X);
}

ScaledVectorGridFunctionAux::~ScaledVectorGridFunctionAux() {
  if (a != nullptr)
    delete a;
  if (a_mixed != nullptr)
    delete a_mixed;
  if (mixed_mat != nullptr)
    delete mixed_mat;
  if (a_mat != nullptr)
    delete a_mat;
  if (solver != nullptr)
    delete solver;
}

} // namespace hephaestus

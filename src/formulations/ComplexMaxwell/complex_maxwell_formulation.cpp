#include "complex_maxwell_formulation.hpp"

namespace hephaestus {

ComplexMaxwellOperator::ComplexMaxwellOperator(
    hephaestus::Problem &problem, const std::string &h_curl_var_complex_name,
    const std::string &h_curl_var_real_name,
    const std::string &h_curl_var_imag_name,
    const std::string &stiffness_coef_name, const std::string &mass_coef_name,
    const std::string &loss_coef_name)
    : ProblemOperator(problem),
      _h_curl_var_complex_name(h_curl_var_complex_name),
      _h_curl_var_real_name(h_curl_var_real_name),
      _h_curl_var_imag_name(h_curl_var_imag_name),
      _stiffness_coef_name(stiffness_coef_name),
      _mass_coef_name(mass_coef_name), _loss_coef_name(loss_coef_name) {}

void ComplexMaxwellOperator::SetGridFunctions() {
  trial_var_names.push_back(_h_curl_var_real_name);
  trial_var_names.push_back(_h_curl_var_imag_name);

  ProblemOperator::SetGridFunctions();

  u_ = new mfem::ParComplexGridFunction(trial_variables.at(0)->ParFESpace());
  *u_ = std::complex(0.0, 0.0);
};

void ComplexMaxwellOperator::Init(mfem::Vector &X) {
  ProblemOperator::Init(X);

  stiffCoef_ = _problem.coefficients.scalars.Get(_stiffness_coef_name);
  massCoef_ = _problem.coefficients.scalars.Get(_mass_coef_name);
  lossCoef_ = _problem.coefficients.scalars.Get(_loss_coef_name);
}

void ComplexMaxwellOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Vector U, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  mfem::ParSesquilinearForm a1_(u_->ParFESpace(), conv_);
  a1_.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_), NULL);
  if (massCoef_) {
    a1_.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_), NULL);
  }
  if (lossCoef_) {
    a1_.AddDomainIntegrator(NULL, new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  mfem::ParLinearForm b1_real_(u_->ParFESpace());
  mfem::ParLinearForm b1_imag_(u_->ParFESpace());
  b1_real_ = 0.0;
  b1_imag_ = 0.0;

  _problem.sources.Apply(&b1_real_);

  mfem::ParComplexLinearForm b1_(u_->ParFESpace(), conv_);
  _problem.bc_map.applyEssentialBCs(_h_curl_var_complex_name, ess_bdr_tdofs_,
                                    *u_, _problem.pmesh.get());
  _problem.bc_map.applyIntegratedBCs(_h_curl_var_complex_name, b1_,
                                     _problem.pmesh.get());
  _problem.bc_map.applyIntegratedBCs(_h_curl_var_complex_name, a1_,
                                     _problem.pmesh.get());

  a1_.Assemble();
  a1_.Finalize();

  b1_.Assemble();
  b1_.real() += b1_real_;
  b1_.imag() += b1_imag_;

  a1_.FormLinearSystem(ess_bdr_tdofs_, *u_, b1_, A1, U, RHS);

  mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver jacobian_solver(MPI_COMM_WORLD);
  jacobian_solver.SetOperator(A_SuperLU);
  jacobian_solver.Mult(RHS, U);
  delete A1C;

  a1_.RecoverFEMSolution(U, b1_, *u_);

  *_problem.gridfunctions.Get(trial_var_names.at(0)) = u_->real();
  *_problem.gridfunctions.Get(trial_var_names.at(1)) = u_->imag();
}

ComplexMaxwellFormulation::ComplexMaxwellFormulation(
    const std::string &alpha_coef_name, const std::string &beta_coef_name,
    const std::string &zeta_coef_name, const std::string &frequency_coef_name,
    const std::string &h_curl_var_complex_name,
    const std::string &h_curl_var_real_name,
    const std::string &h_curl_var_imag_name)
    : FrequencyDomainEMFormulation(), _alpha_coef_name(alpha_coef_name),
      _beta_coef_name(beta_coef_name), _zeta_coef_name(zeta_coef_name),
      _frequency_coef_name(frequency_coef_name),
      _h_curl_var_complex_name(h_curl_var_complex_name),
      _h_curl_var_real_name(h_curl_var_real_name),
      _h_curl_var_imag_name(h_curl_var_imag_name),
      _mass_coef_name(std::string("maxwell_mass")),
      _loss_coef_name(std::string("maxwell_loss")) {}

void ComplexMaxwellFormulation::ConstructOperator() {
  problem->ss_operator = std::make_unique<hephaestus::ComplexMaxwellOperator>(
      *problem, _h_curl_var_complex_name, _h_curl_var_real_name,
      _h_curl_var_imag_name, _alpha_coef_name, _mass_coef_name,
      _loss_coef_name);
  problem->GetOperator()->SetGridFunctions();
}

void ComplexMaxwellFormulation::RegisterGridFunctions() {
  int &myid = GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = GetProblem()->gridfunctions;
  hephaestus::FESpaces &fespaces = GetProblem()->fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_h_curl_var_real_name)) {
    if (myid == 0) {
      MFEM_WARNING(_h_curl_var_real_name
                   << " not found in gridfunctions: building "
                      "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P1"));
    AddGridFunction(_h_curl_var_real_name, std::string("_HCurlFESpace"));
    AddGridFunction(_h_curl_var_imag_name, std::string("_HCurlFESpace"));
  };
}

void ComplexMaxwellFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = GetProblem()->coefficients;

  if (!coefficients.scalars.Has(_frequency_coef_name)) {
    MFEM_ABORT(_frequency_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has("magnetic_permeability")) {
    MFEM_ABORT("Magnetic permeability coefficient not found.");
  }
  if (!coefficients.scalars.Has(_beta_coef_name)) {
    MFEM_ABORT(_beta_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has(_zeta_coef_name)) {
    MFEM_ABORT(_zeta_coef_name + " coefficient not found.");
  }

  freqCoef = dynamic_cast<mfem::ConstantCoefficient *>(
      coefficients.scalars.Get(_frequency_coef_name));

  // define transformed
  coefficients.scalars.Register(
      "_angular_frequency",
      new mfem::ConstantCoefficient(2.0 * M_PI * freqCoef->constant), true);
  coefficients.scalars.Register(
      "_neg_angular_frequency",
      new mfem::ConstantCoefficient(-2.0 * M_PI * freqCoef->constant), true);
  coefficients.scalars.Register(
      "_angular_frequency_sq",
      new mfem::ConstantCoefficient(pow(2.0 * M_PI * freqCoef->constant, 2)),
      true);
  coefficients.scalars.Register(
      "_neg_angular_frequency_sq",
      new mfem::ConstantCoefficient(-pow(2.0 * M_PI * freqCoef->constant, 2)),
      true);

  coefficients.scalars.Register(
      "_inv_angular_frequency",
      new mfem::RatioCoefficient(
          1.0, *coefficients.scalars.Get("_angular_frequency")),
      true);

  coefficients.scalars.Register(
      _mass_coef_name,
      new mfem::TransformedCoefficient(
          coefficients.scalars.Get("_neg_angular_frequency_sq"),
          coefficients.scalars.Get(_zeta_coef_name), prodFunc),
      true);

  coefficients.scalars.Register(
      _loss_coef_name,
      new mfem::TransformedCoefficient(
          coefficients.scalars.Get("_angular_frequency"),
          coefficients.scalars.Get(_beta_coef_name), prodFunc),
      true);

  coefficients.scalars.Register(
      _alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get("magnetic_permeability"),
          fracFunc),
      true);
}

} // namespace hephaestus

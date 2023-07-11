#include "complex_maxwell_formulation.hpp"

namespace hephaestus {

ComplexMaxwellOperator::ComplexMaxwellOperator(
    mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
    hephaestus::GridFunctions &gridfunctions, hephaestus::BCMap &bc_map,
    hephaestus::Coefficients &coefficients, hephaestus::Sources &sources,
    hephaestus::InputParameters &solver_options)
    : FrequencyDomainEquationSystemOperator(pmesh, fespaces, gridfunctions,
                                            bc_map, coefficients, sources,
                                            solver_options),
      h_curl_var_name(solver_options.GetParam<std::string>("HCurlVarName")),
      stiffness_coef_name(
          solver_options.GetParam<std::string>("StiffnessCoefName")),
      mass_coef_name(solver_options.GetParam<std::string>("MassCoefName")),
      loss_coef_name(solver_options.GetParam<std::string>("LossCoefName")) {}

void ComplexMaxwellOperator::SetGridFunctions() {
  state_var_names.push_back(h_curl_var_name + "_real");
  state_var_names.push_back(h_curl_var_name + "_imag");

  FrequencyDomainEquationSystemOperator::SetGridFunctions();

  u_ = new mfem::ParComplexGridFunction(local_test_vars.at(0)->ParFESpace());
  *u_ = std::complex(0.0, 0.0);
};

void ComplexMaxwellOperator::Init(mfem::Vector &X) {
  FrequencyDomainEquationSystemOperator::Init(X);

  stiffCoef_ = _coefficients.scalars.Get(stiffness_coef_name);
  massCoef_ = _coefficients.scalars.Get(mass_coef_name);
  lossCoef_ = _coefficients.scalars.Get(loss_coef_name);
}

void ComplexMaxwellOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Vector U, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  a1_ = new mfem::ParSesquilinearForm(u_->ParFESpace(), conv_);
  a1_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_), NULL);
  if (massCoef_) {
    a1_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_),
                             NULL);
  }
  if (lossCoef_) {
    a1_->AddDomainIntegrator(NULL,
                             new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  // Volume Current Density
  mfem::Vector j(3);
  j = 0.0;
  jrCoef_ = new mfem::VectorConstantCoefficient(j);
  jiCoef_ = new mfem::VectorConstantCoefficient(j);

  mfem::ParLinearForm *b1_real_ = new mfem::ParLinearForm(u_->ParFESpace());
  mfem::ParLinearForm *b1_imag_ = new mfem::ParLinearForm(u_->ParFESpace());
  *b1_real_ = 0.0;
  *b1_imag_ = 0.0;

  _sources.Apply(b1_real_);

  b1_ = new mfem::ParComplexLinearForm(u_->ParFESpace(), conv_);

  _bc_map.applyEssentialBCs(std::string(h_curl_var_name), ess_bdr_tdofs_, *u_,
                            pmesh_);
  _bc_map.applyIntegratedBCs(std::string(h_curl_var_name), *b1_, pmesh_);
  _bc_map.applyIntegratedBCs(std::string(h_curl_var_name), *a1_, pmesh_);

  a1_->Assemble();
  a1_->Finalize();

  b1_->Assemble();
  b1_->real() += *b1_real_;
  b1_->imag() += *b1_imag_;

  a1_->FormLinearSystem(ess_bdr_tdofs_, *u_, *b1_, A1, U, RHS);

  mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(A_SuperLU);
  solver.Mult(RHS, U);
  delete A1C;

  a1_->RecoverFEMSolution(U, *b1_, *u_);

  *_gridfunctions.Get(state_var_names.at(0)) = u_->real();
  *_gridfunctions.Get(state_var_names.at(1)) = u_->imag();
}

ComplexMaxwellFormulation::ComplexMaxwellFormulation()
    : FrequencyDomainFormulation() {
  frequency_coef_name = std::string("frequency");
  h_curl_var_name = std::string("h_curl_var");
  mass_coef_name = std::string("maxwell_mass");
  loss_coef_name = std::string("maxwell_loss");

  zeta_coef_name = std::string("zeta");
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
}

void ComplexMaxwellFormulation::ConstructOperator() {
  hephaestus::InputParameters &solver_options =
      this->GetProblem()->solver_options;
  solver_options.SetParam("HCurlVarName", h_curl_var_name);
  solver_options.SetParam("StiffnessCoefName", alpha_coef_name);
  solver_options.SetParam("MassCoefName", mass_coef_name);
  solver_options.SetParam("LossCoefName", loss_coef_name);
  this->problem->fd_operator =
      std::make_unique<hephaestus::ComplexMaxwellOperator>(
          *(this->problem->pmesh), this->problem->fespaces,
          this->problem->gridfunctions, this->problem->bc_map,
          this->problem->coefficients, this->problem->sources,
          this->problem->solver_options);
  this->problem->fd_operator->SetGridFunctions();
}

void ComplexMaxwellFormulation::RegisterGridFunctions() {
  int &myid = this->GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::FESpaces &fespaces = this->GetProblem()->fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(h_curl_var_name + "_real")) {
    if (myid == 0) {
      MFEM_WARNING(h_curl_var_name << " not found in gridfunctions: building "
                                      "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P1"));
    AddGridFunction(h_curl_var_name + "_real", std::string("_HCurlFESpace"));
    AddGridFunction(h_curl_var_name + "_imag", std::string("_HCurlFESpace"));
  };
}

void ComplexMaxwellFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;

  if (!coefficients.scalars.Has(frequency_coef_name)) {
    MFEM_ABORT(frequency_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has("magnetic_permeability")) {
    MFEM_ABORT("Magnetic permeability coefficient not found.");
  }
  if (!coefficients.scalars.Has(beta_coef_name)) {
    MFEM_ABORT(beta_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has(zeta_coef_name)) {
    MFEM_ABORT(zeta_coef_name + " coefficient not found.");
  }

  freqCoef = dynamic_cast<mfem::ConstantCoefficient *>(
      coefficients.scalars.Get(frequency_coef_name));

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
      mass_coef_name,
      new mfem::TransformedCoefficient(
          coefficients.scalars.Get("_neg_angular_frequency_sq"),
          coefficients.scalars.Get(zeta_coef_name), prodFunc),
      true);

  coefficients.scalars.Register(
      loss_coef_name,
      new mfem::TransformedCoefficient(
          coefficients.scalars.Get("_angular_frequency"),
          coefficients.scalars.Get(beta_coef_name), prodFunc),
      true);

  coefficients.scalars.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get("magnetic_permeability"),
          fracFunc),
      true);
}

} // namespace hephaestus

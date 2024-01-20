#include "complex_maxwell_formulation.hpp"

#include <utility>

namespace hephaestus
{

ComplexMaxwellOperator::ComplexMaxwellOperator(mfem::ParMesh & pmesh,
                                               hephaestus::FESpaces & fespaces,
                                               hephaestus::GridFunctions & gridfunctions,
                                               hephaestus::BCMap & bc_map,
                                               hephaestus::Coefficients & coefficients,
                                               hephaestus::Sources & sources,
                                               hephaestus::InputParameters & solver_options)
  : EquationSystemOperator(
        pmesh, fespaces, gridfunctions, bc_map, coefficients, sources, solver_options),
    h_curl_var_complex_name(solver_options.GetParam<std::string>("HCurlVarComplexName")),
    h_curl_var_real_name(solver_options.GetParam<std::string>("HCurlVarRealName")),
    h_curl_var_imag_name(solver_options.GetParam<std::string>("HCurlVarImagName")),
    stiffness_coef_name(solver_options.GetParam<std::string>("StiffnessCoefName")),
    mass_coef_name(solver_options.GetParam<std::string>("MassCoefName")),
    loss_coef_name(solver_options.GetParam<std::string>("LossCoefName"))
{
}

void
ComplexMaxwellOperator::SetGridFunctions()
{
  state_var_names.push_back(h_curl_var_real_name);
  state_var_names.push_back(h_curl_var_imag_name);

  EquationSystemOperator::SetGridFunctions();

  u_ = new mfem::ParComplexGridFunction(local_test_vars.at(0)->ParFESpace());
  *u_ = std::complex(0.0, 0.0);
};

void
ComplexMaxwellOperator::Init(mfem::Vector & X)
{
  EquationSystemOperator::Init(X);

  stiffCoef_ = _coefficients.scalars.Get(stiffness_coef_name);
  massCoef_ = _coefficients.scalars.Get(mass_coef_name);
  lossCoef_ = _coefficients.scalars.Get(loss_coef_name);
}

void
ComplexMaxwellOperator::Solve(mfem::Vector & X)
{
  mfem::OperatorHandle A1;
  mfem::Vector U, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  mfem::ParSesquilinearForm a1_(u_->ParFESpace(), conv_);
  a1_.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_), nullptr);
  if (massCoef_)
  {
    a1_.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_), nullptr);
  }
  if (lossCoef_)
  {
    a1_.AddDomainIntegrator(nullptr, new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  mfem::ParLinearForm b1_real_(u_->ParFESpace());
  mfem::ParLinearForm b1_imag_(u_->ParFESpace());
  b1_real_ = 0.0;
  b1_imag_ = 0.0;

  _sources.Apply(&b1_real_);

  mfem::ParComplexLinearForm b1_(u_->ParFESpace(), conv_);
  _bc_map.applyEssentialBCs(h_curl_var_complex_name, ess_bdr_tdofs_, *u_, pmesh_);
  _bc_map.applyIntegratedBCs(h_curl_var_complex_name, b1_, pmesh_);
  _bc_map.applyIntegratedBCs(h_curl_var_complex_name, a1_, pmesh_);

  a1_.Assemble();
  a1_.Finalize();

  b1_.Assemble();
  b1_.real() += b1_real_;
  b1_.imag() += b1_imag_;

  a1_.FormLinearSystem(ess_bdr_tdofs_, *u_, b1_, A1, U, RHS);

  auto * A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  auto A1C = std::unique_ptr<mfem::HypreParMatrix>(A1Z->GetSystemMatrix());

  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(A_SuperLU);
  solver.Mult(RHS, U);

  a1_.RecoverFEMSolution(U, b1_, *u_);

  *_gridfunctions.Get(state_var_names.at(0)) = u_->real();
  *_gridfunctions.Get(state_var_names.at(1)) = u_->imag();
}

ComplexMaxwellFormulation::ComplexMaxwellFormulation(std::string alpha_coef_name,
                                                     std::string beta_coef_name,
                                                     std::string zeta_coef_name,
                                                     std::string frequency_coef_name,
                                                     std::string h_curl_var_complex_name,
                                                     std::string h_curl_var_real_name,
                                                     std::string h_curl_var_imag_name)
  : _alpha_coef_name(std::move(alpha_coef_name)),
    _beta_coef_name(std::move(beta_coef_name)),
    _zeta_coef_name(std::move(zeta_coef_name)),
    _frequency_coef_name(std::move(frequency_coef_name)),
    _h_curl_var_complex_name(std::move(h_curl_var_complex_name)),
    _h_curl_var_real_name(std::move(h_curl_var_real_name)),
    _h_curl_var_imag_name(std::move(h_curl_var_imag_name)),
    _mass_coef_name(std::string("maxwell_mass")),
    _loss_coef_name(std::string("maxwell_loss"))
{
}

void
ComplexMaxwellFormulation::ConstructOperator()
{
  hephaestus::InputParameters & solver_options = GetProblem()->solver_options;
  solver_options.SetParam("HCurlVarComplexName", _h_curl_var_complex_name);
  solver_options.SetParam("HCurlVarRealName", _h_curl_var_real_name);
  solver_options.SetParam("HCurlVarImagName", _h_curl_var_imag_name);
  solver_options.SetParam("StiffnessCoefName", _alpha_coef_name);
  solver_options.SetParam("MassCoefName", _mass_coef_name);
  solver_options.SetParam("LossCoefName", _loss_coef_name);
  problem->eq_sys_operator =
      std::make_unique<hephaestus::ComplexMaxwellOperator>(*(problem->pmesh),
                                                           problem->fespaces,
                                                           problem->gridfunctions,
                                                           problem->bc_map,
                                                           problem->coefficients,
                                                           problem->sources,
                                                           problem->solver_options);
  problem->GetOperator()->SetGridFunctions();
}

void
ComplexMaxwellFormulation::RegisterGridFunctions()
{
  int & myid = GetProblem()->myid_;
  hephaestus::GridFunctions & gridfunctions = GetProblem()->gridfunctions;
  hephaestus::FESpaces & fespaces = GetProblem()->fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_h_curl_var_real_name))
  {
    if (myid == 0)
    {
      MFEM_WARNING(_h_curl_var_real_name << " not found in gridfunctions: building "
                                            "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P1"));
    AddGridFunction(_h_curl_var_real_name, std::string("_HCurlFESpace"));
    AddGridFunction(_h_curl_var_imag_name, std::string("_HCurlFESpace"));
  };
}

void
ComplexMaxwellFormulation::RegisterCoefficients()
{
  hephaestus::Coefficients & coefficients = GetProblem()->coefficients;

  if (!coefficients.scalars.Has(_frequency_coef_name))
  {
    MFEM_ABORT(_frequency_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has("magnetic_permeability"))
  {
    MFEM_ABORT("Magnetic permeability coefficient not found.");
  }
  if (!coefficients.scalars.Has(_beta_coef_name))
  {
    MFEM_ABORT(_beta_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has(_zeta_coef_name))
  {
    MFEM_ABORT(_zeta_coef_name + " coefficient not found.");
  }

  freqCoef =
      dynamic_cast<mfem::ConstantCoefficient *>(coefficients.scalars.Get(_frequency_coef_name));

  // define transformed
  coefficients.scalars.Register(
      "_angular_frequency", new mfem::ConstantCoefficient(2.0 * M_PI * freqCoef->constant), true);
  coefficients.scalars.Register("_neg_angular_frequency",
                                new mfem::ConstantCoefficient(-2.0 * M_PI * freqCoef->constant),
                                true);
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
      new mfem::RatioCoefficient(1.0, *coefficients.scalars.Get("_angular_frequency")),
      true);

  coefficients.scalars.Register(
      _mass_coef_name,
      new mfem::TransformedCoefficient(coefficients.scalars.Get("_neg_angular_frequency_sq"),
                                       coefficients.scalars.Get(_zeta_coef_name),
                                       prodFunc),
      true);

  coefficients.scalars.Register(
      _loss_coef_name,
      new mfem::TransformedCoefficient(coefficients.scalars.Get("_angular_frequency"),
                                       coefficients.scalars.Get(_beta_coef_name),
                                       prodFunc),
      true);

  coefficients.scalars.Register(
      _alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get("magnetic_permeability"), fracFunc),
      true);
}

} // namespace hephaestus

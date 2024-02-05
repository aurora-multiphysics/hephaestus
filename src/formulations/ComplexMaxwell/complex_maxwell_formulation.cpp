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
    _h_curl_var_complex_name(solver_options.GetParam<std::string>("HCurlVarComplexName")),
    _h_curl_var_real_name(solver_options.GetParam<std::string>("HCurlVarRealName")),
    _h_curl_var_imag_name(solver_options.GetParam<std::string>("HCurlVarImagName")),
    _stiffness_coef_name(solver_options.GetParam<std::string>("StiffnessCoefName")),
    _mass_coef_name(solver_options.GetParam<std::string>("MassCoefName")),
    _loss_coef_name(solver_options.GetParam<std::string>("LossCoefName"))
{
}

void
ComplexMaxwellOperator::SetGridFunctions()
{
  _state_var_names.push_back(_h_curl_var_real_name);
  _state_var_names.push_back(_h_curl_var_imag_name);

  EquationSystemOperator::SetGridFunctions();

  _u = std::make_unique<mfem::ParComplexGridFunction>(_local_test_vars.at(0)->ParFESpace());
  *_u = std::complex(0.0, 0.0);
};

void
ComplexMaxwellOperator::Init(mfem::Vector & X)
{
  EquationSystemOperator::Init(X);

  _stiff_coef = _coefficients._scalars.GetPtr(_stiffness_coef_name);

  if (_coefficients._scalars.Has(_mass_coef_name))
    _mass_coef = _coefficients._scalars.GetPtr(_mass_coef_name);
  if (_coefficients._scalars.Has(_loss_coef_name))
    _loss_coef = _coefficients._scalars.GetPtr(_loss_coef_name);
}

void
ComplexMaxwellOperator::Solve(mfem::Vector & X)
{
  mfem::OperatorHandle jac;
  mfem::Vector u, rhs;
  mfem::OperatorHandle pc_op;

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient zero_coef(zero_vec);

  mfem::ParSesquilinearForm sqlf(_u->ParFESpace(), _conv);
  sqlf.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*_stiff_coef), nullptr);
  if (_mass_coef)
  {
    sqlf.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*_mass_coef), nullptr);
  }
  if (_loss_coef)
  {
    sqlf.AddDomainIntegrator(nullptr, new mfem::VectorFEMassIntegrator(*_loss_coef));
  }

  mfem::ParLinearForm lf_real(_u->ParFESpace());
  mfem::ParLinearForm lf_imag(_u->ParFESpace());
  lf_real = 0.0;
  lf_imag = 0.0;

  _sources.Apply(&lf_real);

  mfem::ParComplexLinearForm lf(_u->ParFESpace(), _conv);
  _bc_map.ApplyEssentialBCs(_h_curl_var_complex_name, _ess_bdr_tdofs, *_u, _pmesh);
  _bc_map.ApplyIntegratedBCs(_h_curl_var_complex_name, lf, _pmesh);
  _bc_map.ApplyIntegratedBCs(_h_curl_var_complex_name, sqlf, _pmesh);

  sqlf.Assemble();
  sqlf.Finalize();

  lf.Assemble();
  lf.real() += lf_real;
  lf.imag() += lf_imag;

  sqlf.FormLinearSystem(_ess_bdr_tdofs, *_u, lf, jac, u, rhs);

  auto * jac_z = jac.As<mfem::ComplexHypreParMatrix>();
  auto jac_c = std::unique_ptr<mfem::HypreParMatrix>(jac_z->GetSystemMatrix());

  mfem::SuperLURowLocMatrix a_super_lu(*jac_c);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(a_super_lu);
  solver.Mult(rhs, u);

  sqlf.RecoverFEMSolution(u, lf, *_u);

  _gridfunctions.GetRef(_state_var_names.at(0)) = _u->real();
  _gridfunctions.GetRef(_state_var_names.at(1)) = _u->imag();
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
  hephaestus::InputParameters & solver_options = GetProblem()->_solver_options;
  solver_options.SetParam("HCurlVarComplexName", _h_curl_var_complex_name);
  solver_options.SetParam("HCurlVarRealName", _h_curl_var_real_name);
  solver_options.SetParam("HCurlVarImagName", _h_curl_var_imag_name);
  solver_options.SetParam("StiffnessCoefName", _alpha_coef_name);
  solver_options.SetParam("MassCoefName", _mass_coef_name);
  solver_options.SetParam("LossCoefName", _loss_coef_name);
  _problem->_eq_sys_operator =
      std::make_unique<hephaestus::ComplexMaxwellOperator>(*(_problem->_pmesh),
                                                           _problem->_fespaces,
                                                           _problem->_gridfunctions,
                                                           _problem->_bc_map,
                                                           _problem->_coefficients,
                                                           _problem->_sources,
                                                           _problem->_solver_options);
  _problem->GetOperator()->SetGridFunctions();
}

void
ComplexMaxwellFormulation::RegisterGridFunctions()
{
  int & myid = GetProblem()->_myid;
  hephaestus::GridFunctions & gridfunctions = GetProblem()->_gridfunctions;
  hephaestus::FESpaces & fespaces = GetProblem()->_fespaces;

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
  }
}

void
ComplexMaxwellFormulation::RegisterCoefficients()
{
  hephaestus::Coefficients & coefficients = GetProblem()->_coefficients;

  if (!coefficients._scalars.Has(_frequency_coef_name))
  {
    MFEM_ABORT(_frequency_coef_name + " coefficient not found.");
  }
  if (!coefficients._scalars.Has("magnetic_permeability"))
  {
    MFEM_ABORT("Magnetic permeability coefficient not found.");
  }
  if (!coefficients._scalars.Has(_beta_coef_name))
  {
    MFEM_ABORT(_beta_coef_name + " coefficient not found.");
  }
  if (!coefficients._scalars.Has(_zeta_coef_name))
  {
    MFEM_ABORT(_zeta_coef_name + " coefficient not found.");
  }

  _freq_coef = coefficients._scalars.GetPtr<mfem::ConstantCoefficient>(_frequency_coef_name);

  // define transformed
  coefficients._scalars.Register(
      "_angular_frequency",
      std::make_shared<mfem::ConstantCoefficient>(2.0 * M_PI * _freq_coef->constant));
  coefficients._scalars.Register(
      "_neg_angular_frequency",
      std::make_shared<mfem::ConstantCoefficient>(-2.0 * M_PI * _freq_coef->constant));
  coefficients._scalars.Register(
      "_angular_frequency_sq",
      std::make_shared<mfem::ConstantCoefficient>(pow(2.0 * M_PI * _freq_coef->constant, 2)));
  coefficients._scalars.Register(
      "_neg_angular_frequency_sq",
      std::make_shared<mfem::ConstantCoefficient>(-pow(2.0 * M_PI * _freq_coef->constant, 2)));

  coefficients._scalars.Register("_inv_angular_frequency",
                                 std::make_shared<mfem::RatioCoefficient>(
                                     1.0, coefficients._scalars.GetRef("_angular_frequency")));

  coefficients._scalars.Register(_mass_coef_name,
                                 std::make_shared<mfem::TransformedCoefficient>(
                                     coefficients._scalars.GetPtr("_neg_angular_frequency_sq"),
                                     coefficients._scalars.GetPtr(_zeta_coef_name),
                                     prodFunc));

  coefficients._scalars.Register(_loss_coef_name,
                                 std::make_shared<mfem::TransformedCoefficient>(
                                     coefficients._scalars.GetPtr("_angular_frequency"),
                                     coefficients._scalars.GetPtr(_beta_coef_name),
                                     prodFunc));

  coefficients._scalars.Register(
      _alpha_coef_name,
      std::make_shared<mfem::TransformedCoefficient>(
          &_one_coef, coefficients._scalars.GetPtr("magnetic_permeability"), fracFunc));
}

} // namespace hephaestus

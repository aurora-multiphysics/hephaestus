#include "steady_formulation.hpp"

namespace hephaestus {

void HertzOperator::SetVariables() {
  state_var_names.push_back("electric_field_real");
  state_var_names.push_back("electric_field_imag");

  local_test_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _variables, state_var_names);
  e_ = new mfem::ParComplexGridFunction(local_test_vars.at(0)->ParFESpace());
  *e_ = std::complex(0.0, 0.0);

  // Set operator size and block structure
  block_trueOffsets.SetSize(local_test_vars.size() + 1);
  block_trueOffsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    block_trueOffsets[ind + 1] =
        local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  block_trueOffsets.PartialSum();

  true_offsets.SetSize(local_test_vars.size() + 1);
  true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    true_offsets[ind + 1] = local_test_vars.at(ind)->ParFESpace()->GetVSize();
  }
  true_offsets.PartialSum();

  this->height = true_offsets[local_test_vars.size()];
  this->width = true_offsets[local_test_vars.size()];
  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);

  // Populate vector of active auxiliary variables
  active_aux_var_names.resize(0);
  for (auto &aux_var_name : aux_var_names) {
    if (_variables.Has(aux_var_name)) {
      active_aux_var_names.push_back(aux_var_name);
    }
  }
};

void HertzOperator::Init(mfem::Vector &X) {
  // Define material property coefficients
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
  }

  mfem::Vector a1Vec(port_length_vector, 3);

  kc = M_PI / a1Vec.Norml2();
  double omega_ = 2.0 * M_PI * freq_;
  k0 = omega_ * sqrt(epsilon0_ * mu0_);
  k_ = std::complex<double>(0., sqrt(k0 * k0 - kc * kc));

  omegaCoef_ = new mfem::ConstantCoefficient(2.0 * M_PI * freq_);
  negOmegaCoef_ = new mfem::ConstantCoefficient(-2.0 * M_PI * freq_);
  omega2Coef_ = new mfem::ConstantCoefficient(pow(2.0 * M_PI * freq_, 2));
  negOmega2Coef_ = new mfem::ConstantCoefficient(-pow(2.0 * M_PI * freq_, 2));

  // Robin coefficient, but already handled here.
  etaInvCoef_ = new mfem::ConstantCoefficient(-k_.imag() / (mu0_ * omega_));

  // hephaestus::BoundaryCondition waveguide_ports(std::string("robin_1"),
  //                                               mfem::Array<int>({1, 2}));

  // _equation_system->buildEquationSystem(_bc_map, _sources);
  dbcs.SetSize(1);
  dbcs = 0;
  dbcs[0] = 1;

  // abcs.SetSize(1);
  // abcs[0] = 3;
  // abcs = mfem::Array<int>({2, 3});

  ess_bdr_.SetSize(pmesh_->bdr_attributes.Max());
  if (dbcs.Size() == 1 && dbcs[0] == -1) {
    ess_bdr_ = 1;
  } else {
    ess_bdr_ = 0;
    for (int i = 0; i < dbcs.Size(); i++) {
      ess_bdr_[dbcs[i] - 1] = 1;
    }
  }
  e_->ParFESpace()->GetEssentialTrueDofs(ess_bdr_, ess_bdr_tdofs_);

  hephaestus::VectorFunctionDirichletBC *tangential_E_bc =
      dynamic_cast<hephaestus::VectorFunctionDirichletBC *>(
          _bc_map["tangential_E"]);

  erCoef_ = tangential_E_bc->vec_coeff;
  eiCoef_ = tangential_E_bc->vec_coeff_im;

  // e_->ProjectCoefficient(*erCoef_, *eiCoef_);
  e_->ProjectBdrCoefficientTangent(*erCoef_, *eiCoef_, ess_bdr_);

  // Setup various coefficients
  // Create a coefficient describing the dielectric permittivity
  epsCoef_ = new mfem::ConstantCoefficient(epsilon0_);
  muInvCoef_ = new mfem::ConstantCoefficient(1.0 / mu0_);
  sigmaCoef_ = new mfem::ConstantCoefficient(0.0);

  massCoef_ =
      new mfem::TransformedCoefficient(negOmega2Coef_, epsCoef_, prodFunc);
  posMassCoef_ =
      new mfem::TransformedCoefficient(omega2Coef_, epsCoef_, prodFunc);
  if (sigmaCoef_) {
    lossCoef_ =
        new mfem::TransformedCoefficient(omegaCoef_, sigmaCoef_, prodFunc);
  }
}

void HertzOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Vector E, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  a1_ = new mfem::ParSesquilinearForm(e_->ParFESpace(), conv_);
  a1_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*muInvCoef_), NULL);
  a1_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_), NULL);
  if (lossCoef_) {
    a1_->AddDomainIntegrator(NULL,
                             new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  // Volume Current Density
  mfem::Vector j(3);
  j = 0.0;
  jrCoef_ = new mfem::VectorConstantCoefficient(j);
  jiCoef_ = new mfem::VectorConstantCoefficient(j);

  jd_ = new mfem::ParComplexLinearForm(e_->ParFESpace(), conv_);
  jd_->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*jrCoef_),
                           new mfem::VectorFEDomainLFIntegrator(*jiCoef_));

  // _bc_map.applyEssentialBCs(std::string("electric_field"), ess_bdr_tdofs_,
  // *e_,
  //                           pmesh_);
  _bc_map.applyIntegratedBCs(std::string("electric_field"), *jd_, pmesh_);
  hephaestus::VectorRobinBC *robin_bc;
  robin_bc =
      dynamic_cast<hephaestus::VectorRobinBC *>(_bc_map["WaveguidePortIn"]);
  robin_bc->applyBC(*a1_);
  robin_bc =
      dynamic_cast<hephaestus::VectorRobinBC *>(_bc_map["WaveguidePortOut"]);
  robin_bc->applyBC(*a1_);

  a1_->Assemble();
  a1_->Finalize();

  jd_->Assemble();

  a1_->FormLinearSystem(ess_bdr_tdofs_, *e_, *jd_, A1, E, RHS);

  mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(A_SuperLU);
  solver.Mult(RHS, E);
  delete A1C;

  e_->Distribute(E);

  *_variables.Get(state_var_names.at(0)) = e_->real();
  *_variables.Get(state_var_names.at(1)) = e_->imag();
}

HertzFormulation::HertzFormulation() : fd_operator(NULL), oneCoef(1.0){};

hephaestus::HertzOperator *HertzFormulation::CreateFrequencyDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  fd_operator =
      new hephaestus::HertzOperator(pmesh, order, fespaces, variables, bc_map,
                                    domain_properties, sources, solver_options);
  return fd_operator;
};

void HertzFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);
  std::string h_curl_var_name = std::string("electric_field");
  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_curl_var_name)) {
    if (myid == 0) {
      std::cout
          << h_curl_var_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(&pmesh, 1, pmesh.Dimension()), true);
    variables.Register(h_curl_var_name + "_real",
                       new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")),
                       true);
    variables.Register(h_curl_var_name + "_imag",
                       new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")),
                       true);
  };
};

void HertzFormulation::RegisterCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("alpha") == 0) {
    domain_properties.scalar_property_map["alpha"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("alpha")));
  }
  if (domain_properties.scalar_property_map.count("beta") == 0) {
    domain_properties.scalar_property_map["beta"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("beta")));
  }
}

} // namespace hephaestus
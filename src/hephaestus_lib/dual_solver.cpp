// Solves the equations
// ∇⋅s0 = 0
// ∇×(αv) - βu = s0
// dv/dt = -∇×u

// where
// s0 ∈ H(div) source field
// v ∈ H(div)
// u ∈ H(curl)
// p ∈ H1

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (αv, ∇×u') - (βu, u') - (s0, u') - <(αv) × n, u'> = 0
// (dv/dt, v') + (∇×u, v') = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// dv/dt_{n+1} ∈ H(div)
// u_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (αv_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
// <(αv) × n, u'> = 0
// (dv/dt_{n+1}, v') + (∇×u_{n+1}, v') = 0
// using
// v_{n+1} = v_{n} + dt dv/dt_{n+1} = v_{n} - dt ∇×u_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(u_{n+1}, u') = b1(u')
// dv/dt_{n+1} = -∇×u

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') + (αv_{n}, ∇×u') + <(αdt∇×u_{n+1}) × n, u'>

#include "dual_solver.hpp"

namespace hephaestus {

DualFormulation::DualFormulation() : TransientFormulation() {
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
  CreateEquationSystem();
}

hephaestus::TimeDependentEquationSystem *
DualFormulation::CreateEquationSystem() {
  std::vector<std::string> state_var_names;
  state_var_names.resize(2);
  state_var_names.at(0) = "h_curl_var";
  state_var_names.at(1) = "h_div_var";

  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("VariableNames", state_var_names);
  weak_form_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  // need weak curl kernel
  equation_system = new hephaestus::CurlCurlEquationSystem(weak_form_params);
  return equation_system;
}

hephaestus::TimeDomainEquationSystemOperator *
DualFormulation::CreateTimeDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  td_operator =
      new hephaestus::DualOperator(pmesh, order, fespaces, variables, bc_map,
                                   domain_properties, sources, solver_options);
  td_operator->SetEquationSystem(equation_system);

  return td_operator;
};

void DualFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);

  // Register default ParGridFunctions of state variables if not provided
  std::string u_name = equation_system->var_names.at(0);
  if (!variables.Has(u_name)) {
    if (myid == 0) {
      std::cout
          << u_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(
        u_name, new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")), true);
  };

  // Register default ParGridFunctions of state variables if not provided
  std::string v_name = equation_system->var_names.at(1);
  if (!variables.Has(v_name)) {
    if (myid == 0) {
      std::cout
          << v_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HDivFESpace",
        new mfem::common::RT_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(
        v_name, new mfem::ParGridFunction(fespaces.Get("_HDivFESpace")), true);
  };

  RegisterTimeDerivatives(equation_system->var_names, variables);
};

void DualFormulation::RegisterCoefficients(
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

DualOperator::DualOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : TimeDomainEquationSystemOperator(pmesh, order, fespaces, variables,
                                       bc_map, domain_properties, sources,
                                       solver_options),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HDivFESpace_(
          new mfem::common::RT_ParFESpace(&pmesh, order, pmesh.Dimension())),
      a1(NULL), curl(NULL), weakCurl(NULL),
      u_(mfem::ParGridFunction(HCurlFESpace_)),
      v_(mfem::ParGridFunction(HDivFESpace_)),
      du_(mfem::ParGridFunction(HCurlFESpace_)),
      dv_(mfem::ParGridFunction(HDivFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets.SetSize(3);
  true_offsets[0] = 0;
  true_offsets[1] = HCurlFESpace_->GetVSize();
  true_offsets[2] = HDivFESpace_->GetVSize();
  true_offsets.PartialSum();

  this->height = true_offsets[2];
  this->width = true_offsets[2];

  HYPRE_BigInt size_h1 = H1FESpace_->GlobalTrueVSize();
  HYPRE_BigInt size_nd = HCurlFESpace_->GlobalTrueVSize();
  HYPRE_BigInt size_rt = HDivFESpace_->GlobalTrueVSize();
  if (myid_ == 0) {
    std::cout << "Total number of         DOFs: " << size_h1 + size_nd + size_rt
              << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Total number of H1      DOFs: " << size_h1 << std::endl;
    std::cout << "Total number of H(Curl) DOFs: " << size_nd << std::endl;
    std::cout << "Total number of H(Div)  DOFs: " << size_rt << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }
}

void DualOperator::Init(mfem::Vector &X) {
  RegisterVariables();
  _fespaces.Register("_H1FESpace", H1FESpace_, false);
  _fespaces.Register("_HCurlFESpace", HCurlFESpace_, false);
  _fespaces.Register("_HDivFESpace", HDivFESpace_, false);

  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);

  _sources.Init(_variables, _fespaces, _bc_map, _domain_properties);

  this->buildCurl(alphaCoef); // (αv_{n}, ∇×u')
  b1 = new mfem::ParLinearForm(HCurlFESpace_);
  A1 = new mfem::HypreParMatrix;
  X1 = new mfem::Vector;
  B1 = new mfem::Vector;

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  v_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  u_.ProjectCoefficient(Zero_vec);
  v_.ProjectCoefficient(Zero_vec);
}

// /*
// This is the main computational code that computes dX/dt implicitly
// where X is the state vector containing p, u and v.

// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// dv/dt_{n+1} ∈ H(div)
// u_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (αv_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
// <(α∇×u_{n+1}) × n, u'> = 0
// (dv/dt_{n+1}, v') + (∇×u_{n+1}, v') = 0
// using
// v_{n+1} = v_{n} + dt dv/dt_{n+1} = v_{n} - dt ∇×u_{n+1}
// */
void DualOperator::ImplicitSolve(const double dt, const mfem::Vector &X,
                                 mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  dtCoef.constant = dt;

  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  v_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  du_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[0]);
  dv_.MakeRef(HDivFESpace_, dX_dt, true_offsets[1]);

  _domain_properties.SetTime(this->GetTime());

  //////////////////////////////////////////////////////////////////////////////
  // (αv_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u')

  // <(α∇×u_{n+1}) × n, u'> = 0

  // a1(u_{n+1}, u') = b1(u')
  // a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
  // b1(u') = (s0_{n+1}, u') + (αv_{n}, ∇×u') + <(αdt∇×u_{n+1}) × n, u'>

  // (αv_{n}, ∇×u')
  // v_ is a grid function but weakCurl is not parallel assembled so is OK
  weakCurl->MultTranspose(v_, *b1);

  _sources.Apply(b1);

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  _bc_map.applyIntegratedBCs(u_name, *b1, pmesh_);
  if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    delete dtAlphaCoef;
    dtAlphaCoef =
        new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);

    this->buildA1(betaCoef, dtAlphaCoef);
  }
  a1->FormLinearSystem(ess_tdof_list, J_gf, *b1, *A1, *X1, *B1);

  // We only need to create the solver and preconditioner once
  if (a1_solver == NULL) {
    a1_solver = new hephaestus::DefaultHCurlPCGSolver(_solver_options, *A1,
                                                      HCurlFESpace_);
  }
  a1_solver->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, u_);
  du_ = 0.0;

  // Subtract off contribution from source
  _sources.SubtractSources(&u_);

  // dv/dt_{n+1} = -∇×u
  // note curl maps GF to GF
  curl->Mult(u_, dv_);
  dv_ *= -1.0;
}

void DualOperator::buildA1(mfem::Coefficient *Sigma,
                           mfem::Coefficient *DtMuInv) {
  if (a1 != NULL) {
    delete a1;
  }

  // First create and assemble the bilinear form.  For now we assume the mesh
  // isn't moving, the materials are time independent, and dt is constant. So
  // we only need to do this once.

  a1 = new mfem::ParBilinearForm(HCurlFESpace_);
  a1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  a1->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*DtMuInv));
  a1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.

  dt_A1 = dtCoef.constant;
}

void DualOperator::buildCurl(mfem::Coefficient *MuInv) {
  if (curl != NULL) {
    delete curl;
  }
  if (weakCurl != NULL) {
    delete weakCurl;
  }

  curl = new mfem::ParDiscreteLinearOperator(HCurlFESpace_, HDivFESpace_);
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  curl->Assemble();

  weakCurl = new mfem::ParMixedBilinearForm(HCurlFESpace_, HDivFESpace_);
  weakCurl->AddDomainIntegrator(new mfem::VectorFECurlIntegrator(*MuInv));
  weakCurl->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void DualOperator::RegisterVariables() {
  u_name = "h_curl_var";
  u_display_name = "H(Curl) variable";

  v_name = "h_div_var";
  v_display_name = "H(Div) variable";

  _variables.Register(u_name, &u_, false);
  _variables.Register(v_name, &v_, false);
}

void DualOperator::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("alpha") == 0) {
    domain_properties.scalar_property_map["alpha"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("alpha")));
  }
  if (domain_properties.scalar_property_map.count("beta") == 0) {
    domain_properties.scalar_property_map["beta"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("beta")));
  }
  alphaCoef = domain_properties.scalar_property_map["alpha"];
  betaCoef = domain_properties.scalar_property_map["beta"];
}

} // namespace hephaestus

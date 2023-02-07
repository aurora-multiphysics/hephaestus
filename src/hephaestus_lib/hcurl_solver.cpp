// Solves the equations
// ∇⋅s0 = 0
// ∇×(α∇×u) + βdu/dt = s0

// where
// s0 ∈ H(div) source field
// u ∈ H(curl)
// p ∈ H1

// Dirichlet boundaries constrain du/dt
// Integrated boundaries constrain (α∇×u) × n

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (α∇×u, ∇×u') + (βdu/dt, u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// du/dt_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
// - (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
// using
// u_{n+1} = u_{n} + dt du/dt_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(du/dt_{n+1}, u') = b1(u')

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>

// Formulation should add equation system(s) and required variables
// - Add state equation system
// Formulation - factory class for building:
// -- state variables if not present
// -- auxkernels if auxvariables exist
// -- references to required coefficients
// -- TimeDependentOperator for solve

#include "hcurl_solver.hpp"

namespace hephaestus {

HCurlSolver::HCurlSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _order(order),
      _fespaces(fespaces), _variables(variables), _bc_map(bc_map),
      _sources(sources), _domain_properties(domain_properties),
      _solver_options(solver_options), a1_solver(NULL), u_(NULL), du_(NULL) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  state_var_names.resize(1);
  state_var_names.at(0) = "h_curl_var";

  aux_var_names.resize(1);
  aux_var_names.at(0) = "curl h_curl_var";
}

void HCurlSolver::Init(mfem::Vector &X) {
  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);

  _sources.Init(_variables, _fespaces, _bc_map, _domain_properties);

  A1 = new mfem::HypreParMatrix;
  X1 = new mfem::Vector;
  B1 = new mfem::Vector;

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  u_->MakeRef(u_->ParFESpace(), const_cast<mfem::Vector &>(X), true_offsets[0]);
  u_->ProjectCoefficient(Zero_vec);

  _weak_form =
      new hephaestus::CurlCurlWeakForm(u_name, *du_, *u_, alphaCoef, betaCoef);
  _weak_form->buildWeakForm(_bc_map, _sources);
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing p, u and v.

Unknowns
s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
du/dt_{n+1} ∈ H(curl)
p_{n+1} ∈ H1

Fully discretised equations
-(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
(α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
- (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
using
u_{n+1} = u_{n} + dt du/dt_{n+1}
*/
void HCurlSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  u_->MakeRef(u_->ParFESpace(), const_cast<mfem::Vector &>(X), true_offsets[0]);
  du_->MakeRef(u_->ParFESpace(), dX_dt, true_offsets[0]);
  _domain_properties.SetTime(this->GetTime());

  _weak_form->setTimeStep(dt);
  _weak_form->updateWeakForm(_bc_map, _sources);
  _weak_form->FormLinearSystem(*A1, *X1, *B1);
  if (a1_solver == NULL) {
    a1_solver = new hephaestus::DefaultHCurlPCGSolver(_solver_options, *A1,
                                                      _weak_form->test_pfes);
  }
  a1_solver->Mult(*B1, *X1);
  _weak_form->RecoverFEMSolution(*X1, *du_);
}

void HCurlSolver::RegisterVariables() {
  u_name = state_var_names.at(0);

  u_ = _variables.Get(u_name);
  // Register default ParGridFunctions of state variables if not provided
  if (u_ == NULL) {
    if (myid_ == 0) {
      std::cout
          << u_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    _fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(pmesh_, _order, pmesh_->Dimension()),
        true);
    _variables.Register(
        u_name, new mfem::ParGridFunction(_fespaces.Get("_HCurlFESpace")),
        true);
    u_ = _variables.Get(u_name);
  }
  _variables.Register("du", new mfem::ParGridFunction(u_->ParFESpace()), true);
  du_ = _variables.Get("du");

  _fespaces.Register(
      "_H1FESpace",
      new mfem::common::H1_ParFESpace(pmesh_, _order, pmesh_->Dimension()),
      true);

  true_offsets.SetSize(2);
  true_offsets[0] = 0;
  true_offsets[1] = u_->ParFESpace()->GetVSize();
  true_offsets.PartialSum();

  this->height = true_offsets[1];
  this->width = true_offsets[1];

  HYPRE_BigInt size_nd = u_->ParFESpace()->GlobalTrueVSize();
  if (myid_ == 0) {
    std::cout << "Total number of         DOFs: " << size_nd << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Total number of H(Curl) DOFs: " << size_nd << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }

  // Populate vector of active auxiliary variables
  active_aux_var_names.resize(0);
  for (auto &aux_var_name : aux_var_names) {
    if (_variables.Get(aux_var_name) != NULL) {
      active_aux_var_names.push_back(aux_var_name);
    }
  }
}

void HCurlSolver::SetMaterialCoefficients(
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

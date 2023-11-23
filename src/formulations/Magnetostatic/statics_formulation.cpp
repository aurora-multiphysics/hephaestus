// Solves the equations
// ∇⋅s0 = 0
// ∇×(α∇×u) = s0

// where
// s0 ∈ H(div) source field
// u ∈ H(curl)

// Dirichlet boundaries constrain u
// Integrated boundaries constrain (α∇×u) × n

// Weak form (Space discretisation)
// (α∇×u, ∇×u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Time discretisation using implicit scheme:
// Unknowns
// u ∈ H(curl)

// Fully discretised equations
// (α∇×u, ∇×u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Rewritten as
// a1(u, u') = b1(u')

// where
// a1(u, u') = (α∇×u, ∇×u')
// b1(u') = (s0, u') + <(α∇×u) × n, u'>
#include "statics_formulation.hpp"

namespace hephaestus {

StaticsFormulation::StaticsFormulation(const std::string &alpha_coef_name,
                                       const std::string &h_curl_var_name)
    : SteadyStateEMFormulation(), _alpha_coef_name(alpha_coef_name),
      _h_curl_var_name(h_curl_var_name) {}

void StaticsFormulation::ConstructJacobianPreconditioner() {
  std::shared_ptr<mfem::HypreAMS> precond{std::make_shared<mfem::HypreAMS>(
      this->problem->gridfunctions.Get(_h_curl_var_name)->ParFESpace())};
  precond->SetSingularProblem();
  precond->SetPrintLevel(-1);
  this->problem->_jacobian_preconditioner = precond;
}

void StaticsFormulation::ConstructJacobianSolver() {
  std::shared_ptr<mfem::HyprePCG> solver{
      std::make_shared<mfem::HyprePCG>(this->problem->comm)};
  solver->SetTol(1e-16);
  solver->SetMaxIter(1000);
  solver->SetPrintLevel(-1);
  solver->SetPreconditioner(*std::dynamic_pointer_cast<mfem::HypreSolver>(
      this->problem->_jacobian_preconditioner));
  this->problem->_jacobian_solver = solver;
}

void StaticsFormulation::ConstructOperator() {
  this->problem->ss_operator = std::make_unique<hephaestus::StaticsOperator>(
      *(this->problem->pmesh), this->problem->fespaces,
      this->problem->gridfunctions, this->problem->bc_map,
      this->problem->coefficients, this->problem->sources,
      *(this->problem->_jacobian_solver), _h_curl_var_name, _alpha_coef_name);
  this->problem->GetOperator()->SetGridFunctions();
};

void StaticsFormulation::RegisterGridFunctions() {
  int &myid = this->GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::FESpaces &fespaces = this->GetProblem()->fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_h_curl_var_name)) {
    if (myid == 0) {
      MFEM_WARNING(_h_curl_var_name << " not found in gridfunctions: building "
                                       "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(_h_curl_var_name, std::string("_HCurlFESpace"));
  };
};

void StaticsFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;
  if (!coefficients.scalars.Has(_alpha_coef_name)) {
    MFEM_ABORT(_alpha_coef_name + " coefficient not found.");
  }
}

StaticsOperator::StaticsOperator(
    mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
    hephaestus::GridFunctions &gridfunctions, hephaestus::BCMap &bc_map,
    hephaestus::Coefficients &coefficients, hephaestus::Sources &sources,
    mfem::Solver &jacobian_solver, const std::string &h_curl_var_name,
    const std::string &stiffness_coef_name)
    : ProblemOperator(pmesh, fespaces, gridfunctions, bc_map, coefficients,
                      sources, jacobian_solver),
      _h_curl_var_name(h_curl_var_name),
      _stiffness_coef_name(stiffness_coef_name) {}

void StaticsOperator::SetGridFunctions() {
  trial_var_names.push_back(_h_curl_var_name);
  ProblemOperator::SetGridFunctions();
};

void StaticsOperator::Init(mfem::Vector &X) {
  ProblemOperator::Init(X);
  stiffCoef_ = _coefficients.scalars.Get(_stiffness_coef_name);
}

/*
This is the main method that solves for u.

Unknowns
s0 ∈ H(div) divergence-free source field
u ∈ H(curl)

Fully discretised equations
(α∇×u, ∇×u') - (s0, u') - <(α∇×u) × n, u'> = 0
*/
void StaticsOperator::Solve(mfem::Vector &X) {
  mfem::ParGridFunction &a_(*trial_variables.at(0));
  a_ = 0.0;
  mfem::ParLinearForm b1_(a_.ParFESpace());
  b1_ = 0.0;
  mfem::Array<int> ess_bdr_tdofs_;
  _bc_map.applyEssentialBCs(_h_curl_var_name, ess_bdr_tdofs_, a_, pmesh_);
  _bc_map.applyIntegratedBCs(_h_curl_var_name, b1_, pmesh_);
  b1_.Assemble();
  _sources.Apply(&b1_);

  mfem::ParBilinearForm a1_(a_.ParFESpace());
  a1_.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_));
  a1_.Assemble();

  mfem::HypreParMatrix CurlMuInvCurl;
  mfem::Vector &A(trueX.GetBlock(0));
  mfem::Vector &RHS(trueRhs.GetBlock(0));
  a1_.FormLinearSystem(ess_bdr_tdofs_, a_, b1_, CurlMuInvCurl, A, RHS);

  getJacobianSolver()->SetOperator(CurlMuInvCurl);
  getJacobianSolver()->Mult(RHS, A);

  // hephaestus::DefaultHCurlPCGSolver _jacobian_solver(
  //     _solver_options, CurlMuInvCurl, a_.ParFESpace());
  // _jacobian_solver.Mult(RHS, A);
  a1_.RecoverFEMSolution(A, b1_, a_);
}

void StaticsOperator::buildJacobianSolver() {
  // getJacobianSolver()->SetOperator(CurlMuInvCurl);
}

} // namespace hephaestus

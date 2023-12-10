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
      problem->gridfunctions.Get(_h_curl_var_name)->ParFESpace())};
  precond->SetSingularProblem();
  precond->SetPrintLevel(-1);
  problem->jacobian_preconditioner = precond;
}

void StaticsFormulation::ConstructJacobianSolver() {
  std::shared_ptr<mfem::HypreFGMRES> solver{
      std::make_shared<mfem::HypreFGMRES>(problem->comm)};
  solver->SetTol(1e-16);
  solver->SetMaxIter(100);
  solver->SetKDim(10);
  solver->SetPrintLevel(-1);
  solver->SetPreconditioner(*std::dynamic_pointer_cast<mfem::HypreSolver>(
      problem->jacobian_preconditioner));
  problem->jacobian_solver = solver;
}

void StaticsFormulation::ConstructOperator() {
  problem->ss_operator = std::make_unique<hephaestus::StaticsOperator>(
      *problem, _h_curl_var_name, _alpha_coef_name);
  problem->GetOperator()->SetGridFunctions();
};

void StaticsFormulation::RegisterGridFunctions() {
  int &myid = GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = GetProblem()->gridfunctions;
  hephaestus::FESpaces &fespaces = GetProblem()->fespaces;

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
  hephaestus::Coefficients &coefficients = GetProblem()->coefficients;
  if (!coefficients.scalars.Has(_alpha_coef_name)) {
    MFEM_ABORT(_alpha_coef_name + " coefficient not found.");
  }
}

StaticsOperator::StaticsOperator(hephaestus::Problem &problem,
                                 const std::string &h_curl_var_name,
                                 const std::string &stiffness_coef_name)
    : ProblemOperator(problem), _h_curl_var_name(h_curl_var_name),
      _stiffness_coef_name(stiffness_coef_name) {}

void StaticsOperator::SetGridFunctions() {
  trial_var_names.push_back(_h_curl_var_name);
  ProblemOperator::SetGridFunctions();
};

void StaticsOperator::Init(mfem::Vector &X) {
  ProblemOperator::Init(X);
  stiffCoef_ = _problem.coefficients.scalars.Get(_stiffness_coef_name);
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
  _problem.bc_map.applyEssentialBCs(_h_curl_var_name, ess_bdr_tdofs_, a_,
                                    _problem.pmesh.get());
  _problem.bc_map.applyIntegratedBCs(_h_curl_var_name, b1_,
                                     _problem.pmesh.get());
  b1_.Assemble();
  _problem.sources.Apply(&b1_);

  mfem::ParBilinearForm a1_(a_.ParFESpace());
  a1_.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_));
  a1_.Assemble();
  a1_.Finalize();
  mfem::HypreParMatrix CurlMuInvCurl;
  mfem::HypreParVector A(a_.ParFESpace());
  mfem::HypreParVector RHS(a_.ParFESpace());
  a1_.FormLinearSystem(ess_bdr_tdofs_, a_, b1_, CurlMuInvCurl, A, RHS);

  // nonlinear_solver->SetSolver(*jacobian_solver);
  // nonlinear_solver->SetOperator(CurlMuInvCurl);
  // nonlinear_solver->Mult(RHS, A);

  _problem.jacobian_solver->SetOperator(CurlMuInvCurl);
  _problem.jacobian_solver->Mult(RHS, A);
  a1_.RecoverFEMSolution(A, b1_, a_);
}

} // namespace hephaestus

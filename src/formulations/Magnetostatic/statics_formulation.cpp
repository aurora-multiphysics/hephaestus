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

bool isInDomain(const int el, const mfem::Array<int> &dom,
                                  const mfem::ParMesh *mesh) {

  // This is for ghost elements
  if (el < 0)
    return false;

  bool verify = false;

  for (auto sd : dom) {
    if (mesh->GetAttribute(el) == sd)
      verify = true;
  }

  return verify;
}

namespace hephaestus {

StaticsFormulation::StaticsFormulation(const std::string &alpha_coef_name,
                                       const std::string &h_curl_var_name)
    : SteadyStateEMFormulation(), _alpha_coef_name(alpha_coef_name),
      _h_curl_var_name(h_curl_var_name) {}

void StaticsFormulation::ConstructOperator() {
  hephaestus::InputParameters &solver_options = GetProblem()->solver_options;
  solver_options.SetParam("HCurlVarName", _h_curl_var_name);
  solver_options.SetParam("StiffnessCoefName", _alpha_coef_name);
  problem->eq_sys_operator = std::make_unique<hephaestus::StaticsOperator>(
      *(problem->pmesh), problem->fespaces, problem->gridfunctions,
      problem->bc_map, problem->coefficients, problem->sources,
      problem->solver_options);
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

StaticsOperator::StaticsOperator(mfem::ParMesh &pmesh,
                                 hephaestus::FESpaces &fespaces,
                                 hephaestus::GridFunctions &gridfunctions,
                                 hephaestus::BCMap &bc_map,
                                 hephaestus::Coefficients &coefficients,
                                 hephaestus::Sources &sources,
                                 hephaestus::InputParameters &solver_options)
    : EquationSystemOperator(pmesh, fespaces, gridfunctions, bc_map,
                             coefficients, sources, solver_options),
      h_curl_var_name(solver_options.GetParam<std::string>("HCurlVarName")),
      stiffness_coef_name(
          solver_options.GetParam<std::string>("StiffnessCoefName")) {}

void StaticsOperator::SetGridFunctions() {
  state_var_names.push_back(h_curl_var_name);
  EquationSystemOperator::SetGridFunctions();
};

void StaticsOperator::Init(mfem::Vector &X) {
  EquationSystemOperator::Init(X);
  stiffCoef_ = _coefficients.scalars.Get(stiffness_coef_name);
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
  mfem::ParGridFunction &a_(*local_test_vars.at(0));
  a_ = 0.0;
  mfem::ParLinearForm b1_(a_.ParFESpace());
  b1_ = 0.0;
  mfem::Array<int> ess_bdr_tdofs_;
  _bc_map.applyEssentialBCs(h_curl_var_name, ess_bdr_tdofs_, a_, pmesh_);
  _bc_map.applyIntegratedBCs(h_curl_var_name, b1_, pmesh_);
  b1_.Assemble();
  _sources.Apply(&b1_);
  mfem::ParBilinearForm a1_(a_.ParFESpace());
  a1_.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_));
  a1_.Assemble();
  a1_.Finalize();
  mfem::HypreParMatrix CurlMuInvCurl;
  mfem::HypreParVector A(a_.ParFESpace());
  mfem::HypreParVector RHS(a_.ParFESpace());
  a1_.FormLinearSystem(ess_bdr_tdofs_, a_, b1_, CurlMuInvCurl, A, RHS);

  // Define and apply a parallel FGMRES solver for AX=B with the AMS
  // preconditioner from hypre.
  /*
  mfem::HypreParVector int_nodes;
  int_nodes.SetSize(a_.ParFESpace()->GetNV());
  mfem::Array<int> coil_domains({3,4,5,6});

  int_nodes = 0;
  mfem::Table* v_to_el = pmesh_->GetVertexToElementTable();
  v_to_el->Print();
  std::cout << "NUMBER OF ELEMENTS = " << pmesh_->GetNE() << std::endl;
  std::cout << "NUMBER OF VERTICES = " << pmesh_->GetNV() << std::endl;
  std::cout << "NUMBER OF VERTICES FE = " << a_.ParFESpace()->GetNV() << std::endl;

  for (int v = 0; v < pmesh_->GetNV(); ++v)
    for (int e = 0; e < pmesh_->GetNE(); ++e)
      std::cout << *v_to_el(v,e) << std::endl;
  //for (int v = 0; v < int_nodes.Size(); ++v){
    
    for (auto e:v_to_el[v]){
      if (isInDomain(e,coil_domains,pmesh_)){
        int_nodes[v] = 1;
        break;
      }
    }
    
    

  //}
  */
  

  hephaestus::DefaultHCurlFGMRESSolver a1_solver(_solver_options, CurlMuInvCurl,
                                                 a_.ParFESpace());
  a1_solver.Mult(RHS, A);
  a1_.RecoverFEMSolution(A, b1_, a_);
}

} // namespace hephaestus

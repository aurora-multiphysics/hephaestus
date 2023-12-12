#include "problem_builder.hpp"

namespace hephaestus {
void ProblemBuilder::SetMesh(std::shared_ptr<mfem::ParMesh> pmesh) {
  GetProblem()->pmesh = pmesh;
  GetProblem()->comm = pmesh->GetComm();
  MPI_Comm_size(pmesh->GetComm(), &(GetProblem()->num_procs_));
  MPI_Comm_rank(pmesh->GetComm(), &(GetProblem()->myid_));
}

void ProblemBuilder::SetFESpaces(hephaestus::FESpaces &fespaces) {
  GetProblem()->fespaces = fespaces;
}

void ProblemBuilder::SetGridFunctions(
    hephaestus::GridFunctions &gridfunctions) {
  GetProblem()->gridfunctions = gridfunctions;
}

void ProblemBuilder::SetBoundaryConditions(hephaestus::BCMap &bc_map) {
  GetProblem()->bc_map = bc_map;
}

void ProblemBuilder::SetAuxSolvers(hephaestus::AuxSolvers &preprocessors) {
  GetProblem()->preprocessors = preprocessors;
}

void ProblemBuilder::SetPostprocessors(hephaestus::AuxSolvers &postprocessors) {
  GetProblem()->postprocessors = postprocessors;
}

void ProblemBuilder::SetSources(hephaestus::Sources &sources) {
  GetProblem()->sources = sources;
}

void ProblemBuilder::SetOutputs(hephaestus::Outputs &outputs) {
  GetProblem()->outputs = outputs;
}

void ProblemBuilder::SetJacobianPreconditioner(
    std::shared_ptr<mfem::Solver> preconditioner) {
  GetProblem()->jacobian_preconditioner = preconditioner;
}

void ProblemBuilder::SetJacobianSolver(
    std::shared_ptr<mfem::Solver> jacobian_solver) {
  GetProblem()->jacobian_solver = jacobian_solver;
}

void ProblemBuilder::SetCoefficients(hephaestus::Coefficients &coefficients) {
  GetProblem()->coefficients = coefficients;
}

void ProblemBuilder::AddFESpace(std::string fespace_name, std::string fec_name,
                                int vdim, int ordering) {
  if (!GetProblem()->fecs.Has(fec_name)) {
    mfem::FiniteElementCollection *fec =
        mfem::FiniteElementCollection::New(fec_name.c_str());
    GetProblem()->fecs.Register(fec_name, fec, true);
  }

  if (!GetProblem()->fespaces.Has(fespace_name)) {
    mfem::ParMesh *pmesh = GetProblem()->pmesh.get();
    if (pmesh == NULL) {
      MFEM_ABORT("ParMesh not found when trying to add " << fespace_name
                                                         << " to fespaces.");
    }
    mfem::ParFiniteElementSpace *pfes = new mfem::ParFiniteElementSpace(
        GetProblem()->pmesh.get(), GetProblem()->fecs.Get(fec_name), vdim,
        ordering);

    GetProblem()->fespaces.Register(fespace_name, pfes, true);
  }
}

void ProblemBuilder::AddGridFunction(std::string gridfunction_name,
                                     std::string fespace_name) {
  mfem::ParFiniteElementSpace *fespace(
      GetProblem()->fespaces.Get(fespace_name));
  if (fespace == NULL) {
    MFEM_ABORT(
        "FESpace " << fespace_name
                   << " not found in fespaces when trying to add "
                   << gridfunction_name
                   << " associated with it into gridfunctions. Please add "
                   << fespace_name
                   << " to fespaces before adding this gridfunction.");
  }
  mfem::ParGridFunction *gridfunc = new mfem::ParGridFunction(fespace);
  *gridfunc = 0.0;

  GetProblem()->gridfunctions.Register(gridfunction_name, gridfunc, true);
}

void ProblemBuilder::AddBoundaryCondition(std::string bc_name,
                                          hephaestus::BoundaryCondition *bc,
                                          bool own_data) {
  GetProblem()->bc_map.Register(bc_name, bc, own_data);
}

void ProblemBuilder::AddAuxSolver(std::string auxsolver_name,
                                  hephaestus::AuxSolver *aux, bool own_data) {
  GetProblem()->preprocessors.Register(auxsolver_name, aux, own_data);
}

void ProblemBuilder::AddPostprocessor(std::string auxsolver_name,
                                      hephaestus::AuxSolver *aux,
                                      bool own_data) {
  GetProblem()->postprocessors.Register(auxsolver_name, aux, own_data);
}

void ProblemBuilder::AddSource(std::string source_name,
                               hephaestus::Source *source, bool own_data) {
  GetProblem()->sources.Register(source_name, source, own_data);
}

void ProblemBuilder::ConstructJacobianPreconditioner() {
  std::shared_ptr<mfem::HypreBoomerAMG> precond{
      std::make_shared<mfem::HypreBoomerAMG>()};
  precond->SetPrintLevel(-1);
  GetProblem()->jacobian_preconditioner = precond;
}

void ProblemBuilder::ConstructJacobianSolver() {
  std::shared_ptr<mfem::HypreGMRES> solver{
      std::make_shared<mfem::HypreGMRES>(GetProblem()->comm)};
  solver->SetTol(1e-16);
  solver->SetMaxIter(1000);
  solver->SetPrintLevel(-1);
  solver->SetPreconditioner(*std::dynamic_pointer_cast<mfem::HypreSolver>(
      GetProblem()->jacobian_preconditioner));
  GetProblem()->jacobian_solver = solver;
}

void ProblemBuilder::ConstructNonlinearSolver() {
  std::shared_ptr<mfem::NewtonSolver> nl_solver{
      std::make_shared<mfem::NewtonSolver>(GetProblem()->comm)};
  // Defaults to one iteration, without further nonlinear iterations
  nl_solver->SetRelTol(0.0);
  nl_solver->SetAbsTol(0.0);
  nl_solver->SetMaxIter(1);
  GetProblem()->nonlinear_solver = nl_solver;
}

void ProblemBuilder::InitializeAuxSolvers() {
  GetProblem()->preprocessors.Init(GetProblem()->gridfunctions,
                                   GetProblem()->coefficients);
  GetProblem()->postprocessors.Init(GetProblem()->gridfunctions,
                                    GetProblem()->coefficients);
}

void ProblemBuilder::InitializeOutputs() {
  GetProblem()->outputs.Init(GetProblem()->gridfunctions);
}

} // namespace hephaestus

#include "problem_builder.hpp"

namespace hephaestus {

Problem::~Problem() {
  gridfunctions.DeleteData(true);
  fespaces.DeleteData(true);
}

void ProblemBuilder::SetMesh(std::shared_ptr<mfem::ParMesh> pmesh) {
  GetProblem()->pmesh = pmesh;
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

void ProblemBuilder::SetSolverOptions(
    hephaestus::InputParameters &solver_options) {
  GetProblem()->solver_options = solver_options;
}

void ProblemBuilder::SetCoefficients(hephaestus::Coefficients &coefficients) {
  GetProblem()->coefficients = coefficients;
}

void ProblemBuilder::AddFESpace(std::string fespace_name, std::string fec_name,
                                int vdim, int ordering) {
  if (GetProblem()->fespaces.Has(fespace_name)) {
    const std::string error_message =
        "A fespace with the name " + fespace_name +
        " has already been added to the problem fespaces.";
    mfem::mfem_error(error_message.c_str());
  }
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
  if (GetProblem()->gridfunctions.Has(gridfunction_name)) {
    const std::string error_message =
        "A gridfunction with the name " + gridfunction_name +
        " has already been added to the problem gridfunctions.";
    mfem::mfem_error(error_message.c_str());
  }
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
  if (GetProblem()->bc_map.Has(bc_name)) {
    const std::string error_message =
        "A boundary condition with the name " + bc_name +
        " has already been added to the problem boundary conditions.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->bc_map.Register(bc_name, bc, own_data);
}

void ProblemBuilder::AddAuxSolver(std::string auxsolver_name,
                                  hephaestus::AuxSolver *aux, bool own_data) {
  if (GetProblem()->preprocessors.Has(auxsolver_name)) {
    const std::string error_message =
        "An auxsolver with the name " + auxsolver_name +
        " has already been added to the problem preprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->preprocessors.Register(auxsolver_name, aux, own_data);
}

void ProblemBuilder::AddPostprocessor(std::string auxsolver_name,
                                      hephaestus::AuxSolver *aux,
                                      bool own_data) {
  if (GetProblem()->postprocessors.Has(auxsolver_name)) {
    const std::string error_message =
        "An auxsolver with the name " + auxsolver_name +
        " has already been added to the problem postprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->postprocessors.Register(auxsolver_name, aux, own_data);
}

void ProblemBuilder::AddSource(std::string source_name,
                               hephaestus::Source *source, bool own_data) {
  if (GetProblem()->sources.Has(source_name)) {
    const std::string error_message =
        "A source with the name " + source_name +
        " has already been added to the problem sources.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->sources.Register(source_name, source, own_data);
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

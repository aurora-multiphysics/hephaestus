#include "problem_builder.hpp"

namespace hephaestus {
void ProblemBuilder::SetMesh(std::shared_ptr<mfem::ParMesh> pmesh) {
  this->GetProblem()->pmesh = pmesh;
  MPI_Comm_rank(pmesh->GetComm(), &(this->GetProblem()->myid_));
};

void ProblemBuilder::SetFESpaces(hephaestus::FESpaces &fespaces) {
  this->GetProblem()->fespaces = fespaces;
};

void ProblemBuilder::SetGridFunctions(
    hephaestus::GridFunctions &gridfunctions) {
  this->GetProblem()->gridfunctions = gridfunctions;
};

void ProblemBuilder::SetBoundaryConditions(hephaestus::BCMap &bc_map) {
  this->GetProblem()->bc_map = bc_map;
};

void ProblemBuilder::SetAuxSolvers(hephaestus::AuxSolvers &preprocessors) {
  this->GetProblem()->preprocessors = preprocessors;
};

void ProblemBuilder::SetPostprocessors(hephaestus::AuxSolvers &postprocessors) {
  this->GetProblem()->postprocessors = postprocessors;
};

void ProblemBuilder::SetSources(hephaestus::Sources &sources) {
  this->GetProblem()->sources = sources;
};

void ProblemBuilder::SetOutputs(hephaestus::Outputs &outputs) {
  this->GetProblem()->outputs = outputs;
};

void ProblemBuilder::SetSolverOptions(
    hephaestus::InputParameters &solver_options) {
  this->GetProblem()->solver_options = solver_options;
};

void ProblemBuilder::SetCoefficients(
    hephaestus::Coefficients &domain_properties) {
  this->GetProblem()->domain_properties = domain_properties;
};

void ProblemBuilder::AddFESpace(std::string fespace_name, std::string fec_name,
                                int vdim, int ordering) {
  if (!this->GetProblem()->fecs.Has(fec_name)) {
    mfem::FiniteElementCollection *fec =
        mfem::FiniteElementCollection::New(fec_name.c_str());
    this->GetProblem()->fecs.Register(fec_name, fec, true);
  }

  if (!this->GetProblem()->fespaces.Has(fespace_name)) {
    mfem::ParMesh *pmesh = this->GetProblem()->pmesh.get();
    if (pmesh == NULL) {
      MFEM_ABORT("ParMesh not found when trying to add " << fespace_name
                                                         << " to fespaces.");
    }
    mfem::ParFiniteElementSpace *pfes = new mfem::ParFiniteElementSpace(
        this->GetProblem()->pmesh.get(), this->GetProblem()->fecs.Get(fec_name),
        vdim, ordering);

    this->GetProblem()->fespaces.Register(fespace_name, pfes, true);
  }
};

void ProblemBuilder::AddGridFunction(std::string gridfunction_name,
                                     std::string fespace_name) {
  mfem::ParFiniteElementSpace *fespace(
      this->GetProblem()->fespaces.Get(fespace_name));
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

  this->GetProblem()->gridfunctions.Register(gridfunction_name, gridfunc, true);
};

void ProblemBuilder::InitializePostprocessors() {
  this->GetProblem()->postprocessors.Init(
      this->GetProblem()->gridfunctions, this->GetProblem()->domain_properties);
};

} // namespace hephaestus

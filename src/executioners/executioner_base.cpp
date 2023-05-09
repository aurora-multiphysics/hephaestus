#include "executioner.hpp"

namespace hephaestus {

ExecutionerBase::ExecutionerBase(const hephaestus::InputParameters &params)
    : visualization(params.GetOptionalParam<bool>("UseGLVis", false)) {
  // Read in key objects for solve
  pmesh = new mfem::ParMesh(params.GetParam<mfem::ParMesh>("Mesh"));
  bc_map = new hephaestus::BCMap(
      params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  domain_properties = new hephaestus::DomainProperties(
      params.GetParam<hephaestus::DomainProperties>("DomainProperties"));
  fespaces = new hephaestus::FESpaces(
      params.GetParam<hephaestus::FESpaces>("FESpaces"));
  gridfunctions = new hephaestus::GridFunctions(
      params.GetParam<hephaestus::GridFunctions>("GridFunctions"));
  auxkernels = new hephaestus::AuxKernels(
      params.GetParam<hephaestus::AuxKernels>("AuxKernels"));
  postprocessors = new hephaestus::Postprocessors(
      params.GetParam<hephaestus::Postprocessors>("Postprocessors"));
  sources =
      new hephaestus::Sources(params.GetParam<hephaestus::Sources>("Sources"));
  outputs =
      new hephaestus::Outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  data_collections = new std::map<std::string, mfem::DataCollection *>(
      outputs->data_collections);
  solver_options = new hephaestus::InputParameters(
      params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters()));

  MPI_Comm_size(pmesh->GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh->GetComm(), &myid_);
}

} // namespace hephaestus

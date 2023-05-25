#include "executioner.hpp"

namespace hephaestus {

Problem::Problem(const hephaestus::InputParameters &params)
    : pmesh(mfem::ParMesh(params.GetParam<mfem::ParMesh>("Mesh"))),
      bc_map(hephaestus::BCMap(
          params.GetParam<hephaestus::BCMap>("BoundaryConditions"))),
      domain_properties(hephaestus::DomainProperties(
          params.GetParam<hephaestus::DomainProperties>("DomainProperties"))),
      fespaces(hephaestus::FESpaces(
          params.GetParam<hephaestus::FESpaces>("FESpaces"))),
      gridfunctions(hephaestus::GridFunctions(
          params.GetParam<hephaestus::GridFunctions>("GridFunctions"))),
      auxkernels(hephaestus::AuxKernels(
          params.GetParam<hephaestus::AuxKernels>("AuxKernels"))),
      postprocessors(hephaestus::Postprocessors(
          params.GetParam<hephaestus::Postprocessors>("Postprocessors"))),
      sources(
          hephaestus::Sources(params.GetParam<hephaestus::Sources>("Sources"))),
      outputs(
          hephaestus::Outputs(params.GetParam<hephaestus::Outputs>("Outputs"))),
      data_collections(std::map<std::string, mfem::DataCollection *>(
          outputs.data_collections)),
      solver_options(hephaestus::InputParameters(
          params.GetOptionalParam<hephaestus::InputParameters>(
              "SolverOptions", hephaestus::InputParameters()))) {

  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);
}

Executioner::Executioner(const hephaestus::InputParameters &params)
    : visualization(params.GetOptionalParam<bool>("UseGLVis", false)) {}

} // namespace hephaestus

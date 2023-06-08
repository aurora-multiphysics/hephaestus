#include "problem_builder.hpp"

namespace hephaestus {

Problem::Problem(const hephaestus::InputParameters &params)
    : pmesh(std::make_unique<mfem::ParMesh>(
          params.GetParam<mfem::ParMesh>("Mesh"))),
      bc_map(hephaestus::BCMap(
          params.GetParam<hephaestus::BCMap>("BoundaryConditions"))),
      domain_properties(hephaestus::DomainProperties(
          params.GetParam<hephaestus::DomainProperties>("DomainProperties"))),
      fespaces(hephaestus::FESpaces(
          params.GetParam<hephaestus::FESpaces>("FESpaces"))),
      gridfunctions(hephaestus::GridFunctions(
          params.GetParam<hephaestus::GridFunctions>("GridFunctions"))),
      preprocessors(hephaestus::AuxSolvers(
          params.GetParam<hephaestus::AuxSolvers>("PreProcessors"))),
      postprocessors(hephaestus::AuxSolvers(
          params.GetParam<hephaestus::AuxSolvers>("PostProcessors"))),
      sources(
          hephaestus::Sources(params.GetParam<hephaestus::Sources>("Sources"))),
      outputs(
          hephaestus::Outputs(params.GetParam<hephaestus::Outputs>("Outputs"))),
      solver_options(hephaestus::InputParameters(
          params.GetOptionalParam<hephaestus::InputParameters>(
              "SolverOptions", hephaestus::InputParameters()))) {

  MPI_Comm_size(pmesh->GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh->GetComm(), &myid_);
}

} // namespace hephaestus

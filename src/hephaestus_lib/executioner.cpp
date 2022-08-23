#include "executioner.hpp"

namespace hephaestus {

TransientExecutioner::TransientExecutioner(
    const hephaestus::InputParameters &params)
    : dt(params.GetParam<float>("TimeStep")),
      t_initial(params.GetParam<float>("StartTime")),
      t_final(params.GetParam<float>("EndTime")), t(t_initial), vis_steps(1),
      visualization(false), last_step(false) {}

void TransientExecutioner::Solve(
    const hephaestus::InputParameters &params) const {
  // Read in inputs, and initialise solver
  mfem::ParMesh pmesh(params.GetParam<mfem::ParMesh>("Mesh"));
  int order(params.GetParam<int>("Order"));
  hephaestus::BCMap bc_map(
      params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  hephaestus::DomainProperties domain_properties(
      params.GetParam<hephaestus::DomainProperties>("DomainProperties"));
  hephaestus::Variables variables(
      params.GetParam<hephaestus::Variables>("Variables"));
  hephaestus::AuxKernels auxkernels(
      params.GetParam<hephaestus::AuxKernels>("AuxKernels"));
  hephaestus::Postprocessors postprocessors(
      params.GetParam<hephaestus::Postprocessors>("Postprocessors"));
  hephaestus::Outputs outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  std::string formulation_name(params.GetParam<std::string>("FormulationName"));

  hephaestus::TransientFormulation *formulation =
      hephaestus::Factory::createTransientFormulation(
          formulation_name, pmesh, order, variables.gfs, bc_map,
          domain_properties);

  mfem::BlockVector F(formulation->true_offsets); // Vector of dofs
  formulation->Init(F);                           // Set up initial conditions

  formulation->SetTime(t);
  mfem::ODESolver *ode_solver = new mfem::BackwardEulerSolver;
  ode_solver->Init(*formulation);

  variables.Init(pmesh);
  auxkernels.Init(variables.gfs, domain_properties);
  auxkernels.Solve(t);
  postprocessors.Init(variables.gfs, domain_properties);

  // Set up DataCollections to track fields of interest.
  std::map<std::string, mfem::DataCollection *> data_collections(
      outputs.data_collections);
  for (auto const &[name, dc_] : data_collections) {
    formulation->RegisterOutputFields(dc_);
    // Write initial fields to disk
    formulation->WriteOutputFields(dc_, 0);
  }

  // Initialize GLVis visualization and send the initial condition
  // by socket to a GLVis server.
  if (visualization) {
    formulation->InitializeGLVis();
    formulation->DisplayToGLVis();
  }
  // Begin time evolution
  for (int it = 1; !last_step; it++) {
    // Check if current time step is final
    if (t + dt >= t_final - dt / 2) {
      last_step = true;
    }

    // Advance time step.
    ode_solver->Step(F, t, dt);

    // Output data
    if (last_step || (it % vis_steps) == 0) {

      // Output timestep summary to console
      formulation->WriteConsoleSummary(t, it);

      // Make sure all ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh.GetComm());
      auxkernels.Solve(t);

      // Send output fields to GLVis for visualisation
      if (visualization) {
        formulation->DisplayToGLVis();
      }

      // Save output fields at timestep to DataCollections
      for (auto const &[name, dc_] : data_collections) {
        formulation->WriteOutputFields(dc_, it);
      }
      postprocessors.Update(t);
    }
  }
}

} // namespace hephaestus

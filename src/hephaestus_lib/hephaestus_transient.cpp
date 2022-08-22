#include "hephaestus_transient.hpp"

void transient_solve(int argc, char *argv[], hephaestus::Inputs inputs) {
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Read in inputs, and initialise solver
  mfem::ParMesh pmesh = mfem::ParMesh(MPI_COMM_WORLD, inputs.mesh);
  int order = inputs.order;
  hephaestus::BCMap bc_map(inputs.bc_map);
  mfem::NamedFieldsMap<mfem::ParGridFunction> variables;
  hephaestus::DomainProperties domain_properties(inputs.domain_properties);

  hephaestus::TransientFormulation *formulation =
      hephaestus::Factory::createTransientFormulation(inputs.formulation, pmesh,
                                                      order, variables, bc_map,
                                                      domain_properties);

  mfem::BlockVector F(formulation->true_offsets); // Vector of dofs
  formulation->Init(F);                           // Set up initial conditions

  // Set up Executioner
  double t_initial = inputs.executioner.t_initial; // initial time
  double t_final = inputs.executioner.t_final;     // final time
  double dt = inputs.executioner.dt;               // time step
  int vis_steps = 1;
  double t = t_initial; // current time
  formulation->SetTime(t);
  mfem::ODESolver *ode_solver = new mfem::BackwardEulerSolver;
  ode_solver->Init(*formulation);
  std::cout << "\nInitialised solver" << std::endl;

  // Set up DataCollections to track fields of interest.
  std::map<std::string, mfem::DataCollection *> data_collections(
      inputs.outputs.data_collections);
  for (auto const &[name, dc_] : data_collections) {
    formulation->RegisterOutputFields(dc_);
    // Write initial fields to disk
    formulation->WriteOutputFields(dc_, 0);
  }

  // Initialize GLVis visualization and send the initial condition
  // by socket to a GLVis server.
  bool visualization = true;
  if (visualization) {
    formulation->InitializeGLVis();
    formulation->DisplayToGLVis();
  }
  std::cout << "\nStarting timestep" << std::endl;
  // Begin time evolution
  bool last_step = false;
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

      // Send output fields to GLVis for visualisation
      if (visualization) {
        formulation->DisplayToGLVis();
      }

      // Save output fields at timestep to DataCollections
      for (auto const &[name, dc_] : data_collections) {
        formulation->WriteOutputFields(dc_, it);
      }
    }
  }
  if (myid == 0) {
    std::cout << "\nSolved" << std::endl;
  }

  delete ode_solver;
}

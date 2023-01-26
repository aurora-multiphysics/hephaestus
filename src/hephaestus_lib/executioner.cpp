#include "executioner.hpp"

namespace hephaestus {

TransientExecutioner::TransientExecutioner(
    const hephaestus::InputParameters &params)
    : t_step(params.GetParam<float>("TimeStep")),
      t_initial(params.GetParam<float>("StartTime")),
      t_final(params.GetParam<float>("EndTime")), t(t_initial),
      vis_steps(params.GetOptionalParam<int>("VisualisationSteps", 1)),
      visualization(params.GetOptionalParam<bool>("UseGLVis", false)),
      last_step(false) {}

void TransientExecutioner::Init(const hephaestus::InputParameters &params) {
  // Read in inputs, and initialise solver
  pmesh = new mfem::ParMesh(params.GetParam<mfem::ParMesh>("Mesh"));
  int order(params.GetParam<int>("Order"));
  bc_map = new hephaestus::BCMap(
      params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  domain_properties = new hephaestus::DomainProperties(
      params.GetParam<hephaestus::DomainProperties>("DomainProperties"));
  variables = new hephaestus::Variables(
      params.GetParam<hephaestus::Variables>("Variables"));
  auxkernels = new hephaestus::AuxKernels(
      params.GetParam<hephaestus::AuxKernels>("AuxKernels"));
  postprocessors = new hephaestus::Postprocessors(
      params.GetParam<hephaestus::Postprocessors>("Postprocessors"));
  kernels =
      new hephaestus::Kernels(params.GetParam<hephaestus::Kernels>("Kernels"));
  outputs =
      new hephaestus::Outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  data_collections = new std::map<std::string, mfem::DataCollection *>(
      outputs->data_collections);
  solver_options = new hephaestus::InputParameters(
      params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters()));

  std::string formulation_name(params.GetParam<std::string>("FormulationName"));

  formulation = hephaestus::Factory::createTransientFormulation(
      formulation_name, *pmesh, order, variables->fespaces, variables->gfs,
      *bc_map, *domain_properties, *kernels, *solver_options);

  formulation->RegisterVariables();
  variables->Init(*pmesh);
  auxkernels->Init(variables->gfs, *domain_properties);

  F = new mfem::BlockVector(formulation->true_offsets); // Vector of dofs
  formulation->Init(*F); // Set up initial conditions

  formulation->SetTime(t);
  ode_solver = new mfem::BackwardEulerSolver;
  ode_solver->Init(*formulation);

  postprocessors->Init(variables->gfs, *domain_properties);
  auxkernels->Solve(t);

  // Set up DataCollections to track fields of interest.
  for (auto const &[name, dc_] : *data_collections) {
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
}

void TransientExecutioner::Step(double dt, int it) const {
  // Check if current time step is final
  if (t + dt >= t_final - dt / 2) {
    last_step = true;
  }

  // Advance time step.
  ode_solver->Step(*F, t, dt);

  // Output data
  if (last_step || (it % vis_steps) == 0) {

    // Output timestep summary to console
    formulation->WriteConsoleSummary(t, it);

    // Make sure all ranks have sent their 'v' solution before initiating
    // another set of GLVis connections (one from each rank):
    MPI_Barrier(pmesh->GetComm());
    auxkernels->Solve(t);

    // Send output fields to GLVis for visualisation
    if (visualization) {
      formulation->DisplayToGLVis();
    }

    // Save output fields at timestep to DataCollections
    for (auto const &[name, dc_] : *data_collections) {
      formulation->WriteOutputFields(dc_, it);
    }
    postprocessors->Update(t);
  }
}

void TransientExecutioner::Solve() const {
  // Initialise time variables
  t = t_initial;
  last_step = false;

  for (int it = 1; !last_step; it++) {
    Step(t_step, it);
  }
}

} // namespace hephaestus

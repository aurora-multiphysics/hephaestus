//            ---------------------------------------------------
//            ESolver:  Low-Frequency Electrodynamics Simulation
//            ---------------------------------------------------
//
// This miniapp solves low frequency magnetodynamics problems using the E
// formulation.
//
// (ν∇×E, ∇×E') - (σE, E') - (J0, E') - <(ν∇×E) × n, E'> = 0
// -(J0, ∇ V') + <n.J, V'> = 0

#pragma once
#include "../common/pfem_extras.hpp"
#include "e_solver.hpp"
#include "inputs.hpp"

void e_solve(int argc, char *argv[], hephaestus::Inputs inputs) {

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Read the serial mesh from the given mesh file on all processors. We can
  // handle triangular, quadrilateral, tetrahedral and hexahedral meshes
  // with the same code.
  mfem::Mesh *mesh = new mfem::Mesh(inputs.mesh);
  int dim = mesh->Dimension();
  int max_attr = mesh->bdr_attributes.Max();
  mesh->EnsureNCMesh(); // Required for mesh refinement
  mfem::ParMesh pmesh = mfem::ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;

  int order = inputs.order;
  hephaestus::BCMap bc_map(inputs.bc_map);
  hephaestus::DomainProperties domain_properties(inputs.domain_properties);
  hephaestus::ESolver esolver(pmesh, order, bc_map, domain_properties);

  mfem::BlockVector F(esolver.true_offsets);

  esolver.Init(F);

  mfem::ODESolver *ode_solver = new mfem::BackwardEulerSolver;
  bool visualization = true;
  bool visit = true;

  // Initialize GLVis visualization
  if (visualization) {
    esolver.InitializeGLVis();
  }

  // Initialize VisIt visualization
  mfem::VisItDataCollection visit_dc("AV-Parallel", &pmesh);

  if (visit) {
    esolver.RegisterVisItFields(visit_dc);
  }

  // Write initial fields to disk for VisIt
  if (visit) {
    esolver.WriteVisItFields(0);
  }

  // Send the initial condition by socket to a GLVis server.
  if (visualization) {
    esolver.DisplayToGLVis();
  }

  double ti = 0.0; // initial time
  double t_final = inputs.executioner.t_final;
  double dt = inputs.executioner.dt;
  int vis_steps = 1;

  double t = ti;
  esolver.SetTime(t);
  ode_solver->Init(esolver);

  bool last_step = false;
  for (int it = 1; !last_step; it++) {
    if (t + dt >= t_final - dt / 2) {
      last_step = true;
    }

    // F is the vector of dofs, t is the current time, and dt is the time step
    // to advance.
    ode_solver->Step(F, t, dt);

    if (last_step || (it % vis_steps) == 0) {
      double el = esolver.ElectricLosses();
      if (myid == 0) {
        std::cout << std::fixed;
        std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
                  << std::setprecision(3) << t
                  << ",\tdot(E, J) = " << std::setprecision(8) << el
                  << std::endl;
        ;
      }

      // Make sure all ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh.GetComm());

      if (visualization) {
        esolver.DisplayToGLVis();
      }

      if (visit) {
        esolver.WriteVisItFields(it);
      }
    }
  }
  if (myid == 0) {
    std::cout << "\nSolved" << std::endl;
  }
}

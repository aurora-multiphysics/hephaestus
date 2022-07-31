// Based on an H form MMS test provided by Joseph Dean

#include "auxkernels.hpp"
#include "hephaestus_transient.hpp"
#include "postprocessors.hpp"

#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestHFormSource : public testing::Test {
protected:
  static double estimate_convergence_rate(HYPRE_BigInt n_i, HYPRE_BigInt n_imo,
                                          double error_i, double error_imo,
                                          int dim) {
    return std::log(error_i / error_imo) /
           std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
  }

  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }

  static void hdot_bc(const mfem::Vector &x, double t, mfem::Vector &H) {
    H(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI);
    H(1) = 0;
    H(2) = 0;
  }

  static void H_exact_expr(const mfem::Vector &x, double t,
                           mfem::Vector &H_exact) {
    H_exact(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI) * t;
    H_exact(1) = 0;
    H_exact(2) = 0;
  }
  static double sigma_expr(const mfem::Vector &x) {
    // double sigma = 1.0 + 0.25 * cos(M_PI * x(0)) * cos(M_PI * x(1));
    double sigma = 1.0;
    return sigma;
  }
  // Source field
  static void source_field(const mfem::Vector &x, double t, mfem::Vector &f) {
    f(0) = 2 * M_PI * M_PI * t * sin(M_PI * x(1)) * sin(M_PI * x(2)) +
           sigma_expr(x) * sin(M_PI * x(1)) * sin(M_PI * x(2));
    f(1) = 0.0;
    f(2) = 0.0;
  }

  hephaestus::Inputs hform_rod_inputs() {
    hephaestus::Subdomain wire("wire", 1);
    wire.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain air("air", 2);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({wire, air}));

    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *hdotVecCoef =
        new mfem::VectorFunctionCoefficient(3, hdot_bc);
    bc_map["tangential_dHdt"] = new hephaestus::VectorFunctionDirichletBC(
        std::string("magnetic_field"), mfem::Array<int>({1, 2, 3}),
        hdotVecCoef);
    domain_properties.vector_property_map["surface_tangential_dHdt"] =
        hdotVecCoef;
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(1.0);

    mfem::VectorFunctionCoefficient *dBdtSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_field);
    domain_properties.vector_property_map["source"] = dBdtSrcCoef;

    bc_map["ground_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("magnetic_potential"), mfem::Array<int>({1, 2, 3}),
        new mfem::FunctionCoefficient(potential_ground));

    hephaestus::Executioner executioner(std::string("transient"), 0.05, 0.0,
                                        0.05);
    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("HFormVisIt");
    hephaestus::Outputs outputs(data_collections);
    hephaestus::Inputs inputs(mesh, std::string("HForm"), 2, bc_map,
                              domain_properties, executioner, outputs);
    return inputs;
  }
};

TEST_F(TestHFormSource, CheckRun) {
  hephaestus::Inputs inputs(hform_rod_inputs());
  int num_conv_refinements = 3;
  mfem::VectorFunctionCoefficient H_exact(3, H_exact_expr);
  inputs.domain_properties.vector_property_map["h_exact_coeff"] = &H_exact;

  hephaestus::InputParameters hcurlvarparams;
  hcurlvarparams.SetParam("VariableName",
                          std::string("analytic_magnetic_field"));
  hcurlvarparams.SetParam("FESpaceName", std::string("HCurl"));
  hcurlvarparams.SetParam("FESpaceType", std::string("Nedelec"));
  hcurlvarparams.SetParam("order", 2);
  hcurlvarparams.SetParam("components", 3);

  hephaestus::Variables variables;
  variables.AddVariable(hcurlvarparams);

  hephaestus::InputParameters l2errpostprocparams;
  l2errpostprocparams.SetParam("VariableName", std::string("magnetic_field"));
  l2errpostprocparams.SetParam("VectorCoefficientName",
                               std::string("h_exact_coeff"));
  hephaestus::Postprocessors postprocessors;
  postprocessors.Register(
      "L2ErrorPostprocessor",
      new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams), true);

  hephaestus::InputParameters vectorcoeffauxparams;
  vectorcoeffauxparams.SetParam("VariableName",
                                std::string("analytic_magnetic_field"));
  vectorcoeffauxparams.SetParam("VectorCoefficientName",
                                std::string("h_exact_coeff"));

  hephaestus::AuxKernels auxkernels;
  auxkernels.Register(
      "VectorCoefficientAuxKernel",
      new hephaestus::VectorCoefficientAuxKernel(vectorcoeffauxparams), true);

  for (int par_ref_levels = 0; par_ref_levels < num_conv_refinements;
       ++par_ref_levels) {

    {
      // Read in inputs, and initialise solver
      mfem::ParMesh pmesh = mfem::ParMesh(MPI_COMM_WORLD, inputs.mesh);
      int order = inputs.order;
      hephaestus::BCMap bc_map(inputs.bc_map);
      hephaestus::DomainProperties domain_properties(inputs.domain_properties);

      for (int l = 0; l < par_ref_levels; l++) {
        pmesh.UniformRefinement();
      }

      int myid;
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);

      hephaestus::HFormSolver *formulation = new hephaestus::HFormSolver(
          pmesh, order, variables.gfs, bc_map, domain_properties);

      mfem::BlockVector F(formulation->true_offsets); // Vector of dofs
      formulation->Init(F); // Set up initial conditions

      // Set up Executioner
      double t_initial = inputs.executioner.t_initial; // initial time
      double t_final = inputs.executioner.t_final;     // final time
      double dt = inputs.executioner.dt;               // time step
      int vis_steps = 1;
      double t = t_initial; // current time

      formulation->SetTime(t);
      mfem::ODESolver *ode_solver = new mfem::BackwardEulerSolver;
      ode_solver->Init(*formulation);

      variables.Init(pmesh);
      auxkernels.Init(variables.gfs, domain_properties);
      auxkernels.Solve(t);
      postprocessors.Init(variables.gfs, domain_properties);

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
      delete ode_solver;
    }
  }

  hephaestus::L2ErrorVectorPostprocessor l2errpostprocessor =
      *(dynamic_cast<hephaestus::L2ErrorVectorPostprocessor *>(
          postprocessors.Get("L2ErrorPostprocessor")));

  l2errpostprocessor.times.Print();
  l2errpostprocessor.ndofs.Print();
  l2errpostprocessor.l2_errs.Print();

  double r;
  for (std::size_t i = 1; i < l2errpostprocessor.ndofs.Size(); ++i) {
    r = estimate_convergence_rate(
        l2errpostprocessor.ndofs[i], l2errpostprocessor.ndofs[i - 1],
        l2errpostprocessor.l2_errs[i], l2errpostprocessor.l2_errs[i - 1], 3);
    std::cout << r << std::endl;
  }
}

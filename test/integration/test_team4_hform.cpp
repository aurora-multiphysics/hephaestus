// Test based on TEAM Problem 4
// Implements the TEAM Problem 4 benchmark (FELIX brick)
// available from https://www.compumag.org/wp/team/.
// Reference results available from
// A. Kameari, Results for benchmark calculations of problem 4
// (the felix brick), COMPEL, Vol 7. Num 1 (1988).

#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestTEAM4HForm
{
protected:
  static void ExternaldBdt(const mfem::Vector & xv, double t, mfem::Vector & db_dt)
  {
    // External magnetic flux density B = B0*exp(-t/tau)
    double b0(0.1);     // Initial external magnetic field (T)
    double tau(0.0119); // Time constant (s)

    // Calculate rate of change of external flux density dB/dt
    db_dt(0) = 0.0;
    db_dt(1) = 0.0;
    db_dt(2) = -(b0 / tau) * exp(-t / tau);
  }

  static double ExternaldPsidt(const mfem::Vector & xv, double t)
  {
    // Calculate magnetic potential for external source field
    // External magnetic flux density B = B0*exp(-t/tau)
    mfem::Vector dbext_dt(3);
    ExternaldBdt(xv, t, dbext_dt);

    return xv(2) * dbext_dt(2);
  }

  static void BoundarydHdt(const mfem::Vector & x, double t, mfem::Vector & dh_dt)
  {
    // Assumes boundary is sufficiently far from sources
    dh_dt(0) = 0.0;
    dh_dt(1) = 0.0;
    dh_dt(2) = 0.0;
  }

  /// Returns a unique pointer to the constructed problem builder.
  std::unique_ptr<hephaestus::HFormulation> ConfigureProblem()
  {
    // Create Formulation
    auto problem_builder = std::make_unique<hephaestus::HFormulation>(
        "electric_resistivity", "electric_conductivity", "magnetic_permeability", "magnetic_field");
    // Set Mesh
    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team4_symmetrized.g")).c_str(), 1, 1);
    auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

    problem_builder->SetMesh(pmesh);
    problem_builder->AddFESpace("H1", "H1_3D_P1");
    problem_builder->AddFESpace("HCurl", "ND_3D_P1");
    problem_builder->AddFESpace("HDiv", "RT_3D_P0");
    problem_builder->AddFESpace("Scalar_L2", "L2_3D_P0");
    problem_builder->AddFESpace("Vector_L2", "L2_3D_P0", 3);
    problem_builder->AddGridFunction("magnetic_field", "HCurl");

    problem_builder->AddGridFunction("dmagnetic_potential_dt", "H1");
    problem_builder->AddGridFunction("magnetic_flux_density", "HDiv");
    problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");

    problem_builder->AddGridFunction("current_density", "HDiv");
    problem_builder->RegisterCurrentDensityAux("current_density");

    problem_builder->AddGridFunction("electric_field", "HCurl");
    problem_builder->RegisterElectricFieldAux("electric_field");

    hephaestus::Coefficients coefficients = DefineCoefficients();
    problem_builder->SetCoefficients(coefficients);

    hephaestus::Sources sources = DefineSources();
    problem_builder->SetSources(sources);

    hephaestus::Outputs outputs = DefineOutputs();
    problem_builder->SetOutputs(outputs);

    problem_builder->AddBoundaryCondition(
        "tangential_dhdt_bc",
        std::make_shared<hephaestus::VectorDirichletBC>(
            "dmagnetic_field_dt",
            mfem::Array<int>({1, 2, 5, 6}),
            coefficients._vectors.Get("surface_tangential_dHdt")));
    problem_builder->AddBoundaryCondition(
        "magnetic_potential_bc",
        std::make_shared<hephaestus::ScalarDirichletBC>(
            "dmagnetic_potential_dt",
            mfem::Array<int>({1, 2}),
            coefficients._scalars.Get("magnetic_potential_time_derivative")));

    auto fluxmonitor = std::make_shared<hephaestus::FluxMonitorAux>("current_density", 3);
    fluxmonitor->SetPriority(2);
    problem_builder->AddPostprocessor("FluxMonitor", fluxmonitor);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("AbsTolerance", float(1.0e-20));
    solver_options.SetParam("Tolerance", float(1.0e-20));
    solver_options.SetParam("MaxIter", (unsigned int)500);
    problem_builder->SetSolverOptions(solver_options);

    problem_builder->FinalizeProblem();

    return problem_builder;
  }

  hephaestus::InputParameters DefineExecutionerParameters(hephaestus::TimeDomainProblem & problem)
  {
    hephaestus::InputParameters exec_params;

    exec_params.SetParam("TimeStep", float(0.001));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(0.015));
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("Problem", &problem);

    return exec_params;
  }

  hephaestus::Coefficients DefineCoefficients()
  {
    hephaestus::Subdomain brick("brick", 1);
    brick._scalar_coefficients.Register("electric_conductivity",
                                        std::make_shared<mfem::ConstantCoefficient>(2.538e7));
    hephaestus::Subdomain vacuum("vacuum", 2);
    vacuum._scalar_coefficients.Register("electric_conductivity",
                                         std::make_shared<mfem::ConstantCoefficient>(1.0));

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({brick, vacuum}));
    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

    auto dpsi_dt_coef = std::make_shared<mfem::FunctionCoefficient>(ExternaldPsidt);

    // Register to prevent dpsi_dt_coef being destroyed when it goes out of scope.
    coefficients._scalars.Register("magnetic_potential_time_derivative", std::move(dpsi_dt_coef));

    auto dh_dt_vec_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, BoundarydHdt);
    coefficients._vectors.Register("surface_tangential_dHdt", dh_dt_vec_coef);

    return coefficients;
  }

  hephaestus::Sources DefineSources()
  {
    hephaestus::InputParameters source_solver_options;
    source_solver_options.SetParam("Tolerance", float(1.0e-20));
    source_solver_options.SetParam("MaxIter", (unsigned int)2000);

    hephaestus::Sources sources;
    sources.Register("source",
                     std::make_shared<hephaestus::ScalarPotentialSource>("dhext_dt",
                                                                         "dmagnetic_potential_dt",
                                                                         "HCurl",
                                                                         "H1",
                                                                         "_one",
                                                                         -1.0,
                                                                         source_solver_options));

    return sources;
  }

  hephaestus::Outputs DefineOutputs()
  {
    hephaestus::Outputs outputs;
    outputs.Register("ParaViewDataCollection",
                     std::make_shared<mfem::ParaViewDataCollection>("Team4ParaView"));
    return outputs;
  }

  /// Sets the peak current and time using the flux monitor.
  void ExtractPeakCurrentAndTime(hephaestus::TimeDomainProblem & problem,
                                 double & peak_current,
                                 double & peak_current_time)
  {
    auto flux_monitor = problem._postprocessors.Get<hephaestus::FluxMonitorAux>("FluxMonitor");

    const mfem::real_t min_flux = flux_monitor->_fluxes.Min();

    peak_current = -2.0 * min_flux;
    peak_current_time = 0.0;

    for (int i = 0; i < flux_monitor->_times.Size(); i++)
    {
      const mfem::real_t flux = flux_monitor->_fluxes[i];
      const mfem::real_t flux_time = flux_monitor->_times[i];

      if (flux == min_flux)
      {
        peak_current_time = flux_time;
      }
    }
  }

  /// Checks peak current and time are within expected range.
  void VerifyPeakCurrentAndTimeWithinTolerances(double & peak_current, double & peak_current_time)
  {
    const double min_peak_current = 3.2e3;
    const double max_peak_current = 3.6e3;

    const double min_peak_current_time = 0.010;
    const double max_peak_current_time = 0.012;

    REQUIRE(peak_current > min_peak_current);
    REQUIRE(peak_current < max_peak_current);

    REQUIRE(peak_current_time > min_peak_current_time);
    REQUIRE(peak_current_time < max_peak_current_time);
  }

  /// Calls Reset method of FluxMonitorAux.
  void ResetFluxMonitor(const hephaestus::TimeDomainProblem & problem)
  {
    auto fluxmonitor = problem._postprocessors.Get<hephaestus::FluxMonitorAux>("FluxMonitor");

    fluxmonitor->Reset();
  }
};

TEST_CASE_METHOD(TestTEAM4HForm, "TestTEAM4HForm", "[CheckRun]")
{
  auto problem_builder = ConfigureProblem();
  auto problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params = DefineExecutionerParameters(*problem);

  auto executioner = std::make_unique<hephaestus::TransientExecutioner>(exec_params);
  executioner->Execute();

  double peak_current, peak_current_time;
  ExtractPeakCurrentAndTime(*problem, peak_current, peak_current_time);

  VerifyPeakCurrentAndTimeWithinTolerances(peak_current, peak_current_time);
}

/// Test case for checking "Updates" work as expected following a mesh refinement.
TEST_CASE_METHOD(TestTEAM4HForm, "TestTEAM4HFormMeshUpdates", "[CheckRun][!benchmark]")
{
  auto problem_builder = ConfigureProblem();
  auto problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params = DefineExecutionerParameters(*problem);
  hephaestus::logger.set_level(spdlog::level::info);

  const int imax_refinement = 2;
  for (int irefinement = 0; irefinement < imax_refinement; irefinement++)
  {
    // Refine and reset flux monitor on subsequent runs.
    if (irefinement != 0)
    {
      problem->_pmesh->UniformRefinement();
      problem->Update();

      ResetFluxMonitor(*problem);
    }

    // Create a new executioner.
    hephaestus::TransientExecutioner executioner(exec_params);
    executioner.Execute();

    double peak_current, peak_current_time;
    ExtractPeakCurrentAndTime(*problem, peak_current, peak_current_time);

    hephaestus::logger.info("Refinement Level: {}, Peak current: {}, Time: {} seconds",
                            irefinement,
                            peak_current,
                            peak_current_time);

    VerifyPeakCurrentAndTimeWithinTolerances(peak_current, peak_current_time);
  }
}

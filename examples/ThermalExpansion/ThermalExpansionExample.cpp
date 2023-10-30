#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

hephaestus::BCMap defineBCs(){

  hephaestus::BCMap boundaries;

  boundaries.Register("temperature", new hephaestus::EssentialBC("thermal_boundary_one", {1}, ), );
}


hephaestus::Coefficients defineCoefficients(){
    
  hephaestus::Coefficients coefficients;

  coefficients.scalars.Register("thermal_expansion_coef",
                                new mfem::ConstantCoefficient(0.02),
                                true);

  coefficients.scalars.Register("thermal_conductivity",
                                new mfem::ConstantCoefficient(300),
                                true);

  coefficients.scalars.Register("lame_param",
                                new mfem::ConstantCoefficient(0.02),
                                true);

  coefficients.scalars.Register("shear_modulus",
                                new mfem::ConstantCoefficient(0.02),
                                true);                        

  coefficients.scalars.Register("stress_free_temp",
                                new mfem::ConstantCoefficient(0.02),
                                true);                                                      

  return coefficients;
}

hephaestus::Outputs defineOutputs() {
  std::map<std::string, mfem::DataCollection *> data_collections;
  data_collections["ParaViewDataCollection"] =
      new mfem::ParaViewDataCollection("ThermalExpansionExample");
  hephaestus::Outputs outputs(data_collections);
  return outputs;
}

int main(int argc, char *argv[]) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);

  // Create Formulation
  hephaestus::SteadyStateProblemBuilder *problem_builder = new hephaestus::ThermalExpansionFormulation();
  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./simple_cube.g")).c_str(), 1,
                  1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("H1_3"), std::string("H1_3D_P1"), 3, mfem::Ordering::byVDIM);
  
  problem_builder->AddGridFunction(std::string("temperature"),
                                   std::string("H1"));
  problem_builder->AddGridFunction(std::string("displacement"),
                                   std::string("H1_3"));

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);


  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-5));
  solver_options.SetParam("MaxIter", (unsigned int)1000);
  solver_options.SetParam("PrintLevel", 0);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder);
  sequencer.ConstructOperatorProblem();
  std::unique_ptr<hephaestus::SteadyStateProblem> problem =
      problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", problem.get());
  hephaestus::SteadyStateExecutioner *executioner =
      new hephaestus::SteadyStateExecutioner(exec_params);

  mfem::out << "Created executioner";
  executioner->Init();
  executioner->Execute();

  MPI_Finalize();
}

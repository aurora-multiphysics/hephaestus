#include "problem_builder.hpp"

namespace hephaestus
{

void
Problem::Update()
{
  // 1. Update the fespaces.
  _fespaces.Update();

  // 2. Update the gridfunctions.
  _gridfunctions.Update();

  // 3. Update sources.
  _sources.Update();

  // 4. Update the preprocessors and postprocessors.
  _preprocessors.Update();
  _postprocessors.Update();

  // 5. Rebuild operator (and equation system).
  GetOperator()->Update();
}

void
ProblemBuilder::SetMesh(std::shared_ptr<mfem::ParMesh> pmesh)
{
  logger.info("Setting Mesh");
  GetProblem()->_pmesh = pmesh;
  GetProblem()->_comm = pmesh->GetComm();
  MPI_Comm_size(pmesh->GetComm(), &(GetProblem()->_num_procs));
  MPI_Comm_rank(pmesh->GetComm(), &(GetProblem()->_myid));
}

void
ProblemBuilder::SetFESpaces(hephaestus::FESpaces & fespaces)
{
  logger.info("Setting FE Spaces");
  GetProblem()->_fespaces = fespaces;
}

void
ProblemBuilder::SetGridFunctions(hephaestus::GridFunctions & gridfunctions)
{
  logger.info("Setting GridFunctions");
  GetProblem()->_gridfunctions = gridfunctions;
}

void
ProblemBuilder::SetBoundaryConditions(hephaestus::BCMap & bc_map)
{
  logger.info("Setting Boundary Conditions");
  GetProblem()->_bc_map = bc_map;
}

void
ProblemBuilder::SetAuxSolvers(hephaestus::AuxSolvers & preprocessors)
{
  logger.info("Setting AuxSolvers");
  GetProblem()->_preprocessors = preprocessors;
}

void
ProblemBuilder::SetPostprocessors(hephaestus::AuxSolvers & postprocessors)
{
  logger.info("Setting Postprocessors");
  GetProblem()->_postprocessors = postprocessors;
}

void
ProblemBuilder::SetSources(hephaestus::Sources & sources)
{
  logger.info("Setting Sources");
  GetProblem()->_sources = sources;
}

void
ProblemBuilder::SetOutputs(hephaestus::Outputs & outputs)
{
  logger.info("Setting Outputs");
  GetProblem()->_outputs = outputs;
}

void
ProblemBuilder::SetSolverOptions(hephaestus::InputParameters & solver_options)
{
  logger.info("Setting Solver Options");
  GetProblem()->_solver_options = solver_options;
}

void
ProblemBuilder::SetJacobianPreconditioner(std::unique_ptr<mfem::Solver> preconditioner)
{
  GetProblem()->GetOperator()->SetJacobianPreconditioner(std::move(preconditioner));
}

void
ProblemBuilder::SetJacobianSolver(std::unique_ptr<mfem::Solver> jacobian_solver)
{
  GetProblem()->GetOperator()->SetJacobianSolver(std::move(jacobian_solver));
}

void
ProblemBuilder::SetCoefficients(hephaestus::Coefficients & coefficients)
{
  logger.info("Setting Coefficients");
  GetProblem()->_coefficients = coefficients;
}

void
ProblemBuilder::AddFESpace(std::string fespace_name, std::string fec_name, int vdim, int ordering)
{
  logger.info("Adding {} FE Space to problem", fespace_name);
  if (GetProblem()->_fespaces.Has(fespace_name))
  {
    const std::string error_message = "A fespace with the name " + fespace_name +
                                      " has already been added to the problem fespaces.";
    mfem::mfem_error(error_message.c_str());
  }
  if (!GetProblem()->_fecs.Has(fec_name))
  {
    auto fec = std::shared_ptr<mfem::FiniteElementCollection>(
        mfem::FiniteElementCollection::New(fec_name.c_str()));
    GetProblem()->_fecs.Register(fec_name, fec);
  }

  if (!GetProblem()->_fespaces.Has(fespace_name))
  {
    mfem::ParMesh * pmesh = GetProblem()->_pmesh.get();
    if (pmesh == nullptr)
    {
      MFEM_ABORT("ParMesh not found when trying to add " << fespace_name << " to fespaces.");
    }
    auto pfes = std::make_shared<mfem::ParFiniteElementSpace>(
        GetProblem()->_pmesh.get(), GetProblem()->_fecs.Get(fec_name), vdim, ordering);

    GetProblem()->_fespaces.Register(fespace_name, std::move(pfes));
  }
}

void
ProblemBuilder::AddGridFunction(std::string gridfunction_name, std::string fespace_name)
{
  if (GetProblem()->_gridfunctions.Has(gridfunction_name))
  {
    const std::string error_message = "A gridfunction with the name " + gridfunction_name +
                                      " has already been added to the problem gridfunctions.";
    mfem::mfem_error(error_message.c_str());
  }

  if (!GetProblem()->_fespaces.Has(fespace_name))
  {
    MFEM_ABORT("FESpace " << fespace_name << " not found in fespaces when trying to add "
                          << gridfunction_name
                          << " associated with it into gridfunctions. Please add " << fespace_name
                          << " to fespaces before adding this gridfunction.");
  }

  auto fespace = GetProblem()->_fespaces.Get(fespace_name);

  auto gridfunc = std::make_shared<mfem::ParGridFunction>(fespace);
  *gridfunc = 0.0;

  GetProblem()->_gridfunctions.Register(gridfunction_name, std::move(gridfunc));
}

void
ProblemBuilder::AddBoundaryCondition(std::string bc_name,
                                     std::shared_ptr<hephaestus::BoundaryCondition> bc)
{
  if (GetProblem()->_bc_map.Has(bc_name))
  {
    const std::string error_message = "A boundary condition with the name " + bc_name +
                                      " has already been added to the problem boundary conditions.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_bc_map.Register(bc_name, std::move(bc));
}

void
ProblemBuilder::AddAuxSolver(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux)
{
  if (GetProblem()->_preprocessors.Has(auxsolver_name))
  {
    const std::string error_message = "An auxsolver with the name " + auxsolver_name +
                                      " has already been added to the problem preprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_preprocessors.Register(auxsolver_name, std::move(aux));
}

void
ProblemBuilder::AddPostprocessor(std::string auxsolver_name,
                                 std::shared_ptr<hephaestus::AuxSolver> aux)
{
  if (GetProblem()->_postprocessors.Has(auxsolver_name))
  {
    const std::string error_message = "An auxsolver with the name " + auxsolver_name +
                                      " has already been added to the problem postprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_postprocessors.Register(auxsolver_name, std::move(aux));
}

void
ProblemBuilder::AddSource(std::string source_name, std::shared_ptr<hephaestus::Source> source)
{
  if (GetProblem()->_sources.Has(source_name))
  {
    const std::string error_message =
        "A source with the name " + source_name + " has already been added to the problem sources.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_sources.Register(source_name, std::move(source));
}

void
ProblemBuilder::ConstructJacobianSolver()
{
  GetProblem()->GetOperator()->ConstructJacobianSolver();
}

void
ProblemBuilder::ConstructNonlinearSolver()
{
  auto nl_solver = std::make_unique<mfem::NewtonSolver>(GetProblem()->_comm);

  // Defaults to one iteration, without further nonlinear iterations
  nl_solver->SetRelTol(0.0);
  nl_solver->SetAbsTol(0.0);
  nl_solver->SetMaxIter(1);

  GetProblem()->GetOperator()->SetNonlinearSolver(std::move(nl_solver));
}

void
ProblemBuilder::InitializeSources()
{
  GetProblem()->_sources.Init(GetProblem()->_gridfunctions,
                              GetProblem()->_fespaces,
                              GetProblem()->_bc_map,
                              GetProblem()->_coefficients);
}

void
ProblemBuilder::InitializeAuxSolvers()
{
  GetProblem()->_preprocessors.Init(GetProblem()->_gridfunctions, GetProblem()->_coefficients);
  GetProblem()->_postprocessors.Init(GetProblem()->_gridfunctions, GetProblem()->_coefficients);
}

void
ProblemBuilder::InitializeOutputs()
{
  GetProblem()->_outputs.Init(GetProblem()->_gridfunctions);
}

void
ProblemBuilder::InitializeOperator()
{
  // Setup initial conditions.
  GetProblem()->GetOperator()->Init();
}

void
ProblemBuilder::FinalizeProblem(bool build_operator)
{
  RegisterFESpaces();
  RegisterGridFunctions();
  RegisterAuxSolvers();
  RegisterCoefficients();

  if (build_operator)
  {
    ConstructOperator();
  }

  InitializeAuxSolvers();
  InitializeSources();
  InitializeOperator();

  ConstructJacobianSolver();
  ConstructNonlinearSolver();

  ConstructTimestepper();
  InitializeOutputs();
}

} // namespace hephaestus

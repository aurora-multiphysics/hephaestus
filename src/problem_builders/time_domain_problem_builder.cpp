#include "time_domain_problem_builder.hpp"

namespace hephaestus {

std::vector<mfem::ParGridFunction *>
TimeDomainProblemBuilder::RegisterTimeDerivatives(
    std::vector<std::string> gridfunction_names,
    hephaestus::GridFunctions &gridfunctions) {
  std::vector<mfem::ParGridFunction *> time_derivatives;

  for (auto &gridfunction_name : gridfunction_names) {
    gridfunctions.Register(
        GetTimeDerivativeName(gridfunction_name),
        new mfem::ParGridFunction(
            gridfunctions.Get(gridfunction_name)->ParFESpace()),
        true);
    time_derivatives.push_back(
        gridfunctions.Get(GetTimeDerivativeName(gridfunction_name)));
  }
  return time_derivatives;
}

void TimeDomainProblemBuilder::RegisterGridFunctions() {
  std::vector<std::string> gridfunction_names;
  for (auto const &[name, gf] : this->problem->gridfunctions) {
    gridfunction_names.push_back(name);
  }
  this->RegisterTimeDerivatives(gridfunction_names,
                                this->problem->gridfunctions);
}

void TimeDomainProblemBuilder::ConstructEquationSystem() {
  hephaestus::InputParameters params;
  this->problem->td_equation_system =
      std::make_unique<hephaestus::TimeDependentEquationSystem>(params);
}

void TimeDomainProblemBuilder::InitializeKernels() {
  this->problem->td_equation_system->Init(
      this->problem->gridfunctions, this->problem->fespaces,
      this->problem->bc_map, this->problem->coefficients);

  this->problem->preprocessors.Init(this->problem->gridfunctions,
                                    this->problem->coefficients);
  this->problem->sources.Init(this->problem->gridfunctions,
                              this->problem->fespaces, this->problem->bc_map,
                              this->problem->coefficients);
}

void TimeDomainProblemBuilder::ConstructJacobianPreconditioner() {
  std::shared_ptr<mfem::HypreBoomerAMG> precond{
      std::make_shared<mfem::HypreBoomerAMG>()};
  precond->SetPrintLevel(-1);
  this->problem->_jacobian_preconditioner = precond;
}

void TimeDomainProblemBuilder::ConstructJacobianSolver() {
  std::shared_ptr<mfem::HypreGMRES> solver{
      std::make_shared<mfem::HypreGMRES>(this->problem->comm)};
  solver->SetTol(1e-16);
  solver->SetMaxIter(1000);
  solver->SetPrintLevel(-1);
  solver->SetPreconditioner(*std::dynamic_pointer_cast<mfem::HypreSolver>(
      this->problem->_jacobian_preconditioner));
  this->problem->_jacobian_solver = solver;
}

void TimeDomainProblemBuilder::ConstructNonlinearSolver() {
  std::shared_ptr<mfem::NewtonSolver> nl_solver{
      std::make_shared<mfem::NewtonSolver>(this->problem->comm)};
  // Defaults to one iteration, without further nonlinear iterations
  nl_solver->SetRelTol(0.0);
  nl_solver->SetAbsTol(0.0);
  nl_solver->SetMaxIter(1);
  this->problem->_nonlinear_solver = nl_solver;
}

void TimeDomainProblemBuilder::ConstructOperator() {
  this->problem->td_operator =
      std::make_unique<hephaestus::TimeDomainProblemOperator>(
          *(this->problem->pmesh), this->problem->fespaces,
          this->problem->gridfunctions, this->problem->bc_map,
          this->problem->coefficients, this->problem->sources,
          *(this->problem->_jacobian_solver),
          *(this->problem->_nonlinear_solver));
  this->problem->td_operator->SetEquationSystem(
      this->problem->td_equation_system.get());
  this->problem->td_operator->SetGridFunctions();
}

void TimeDomainProblemBuilder::ConstructState() {
  this->problem->F = new mfem::BlockVector(
      this->problem->td_operator->true_offsets); // Vector of dofs
  this->problem->td_operator->Init(
      *(this->problem->F)); // Set up initial conditions
  this->problem->td_operator->SetTime(0.0);
}

void TimeDomainProblemBuilder::ConstructTimestepper() {
  this->problem->ode_solver = new mfem::BackwardEulerSolver;
  this->problem->ode_solver->Init(*(this->problem->td_operator));
}

} // namespace hephaestus

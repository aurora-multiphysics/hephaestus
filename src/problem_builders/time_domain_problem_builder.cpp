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

void TimeDomainProblemBuilder::ConstructOperator() {
  this->problem->td_operator =
      std::make_unique<hephaestus::TimeDomainEquationSystemOperator>(
          *(this->problem->pmesh), this->problem->fespaces,
          this->problem->gridfunctions, this->problem->bc_map,
          this->problem->coefficients, this->problem->sources,
          this->problem->solver_options);
  this->problem->td_operator->SetEquationSystem(
      this->problem->td_equation_system.get());
  this->problem->td_operator->SetGridFunctions();
}

void TimeDomainProblemBuilder::ConstructState() {
  this->problem->F = new mfem::BlockVector(
      this->problem->td_operator->true_offsets); // Vector of dofs
  *(problem->F) = 0.0;                           // give initial value
  this->problem->td_operator->Init(
      *(this->problem->F)); // Set up initial conditions
  this->problem->td_operator->SetTime(0.0);
}

void TimeDomainProblemBuilder::ConstructSolver() {
  this->problem->ode_solver = new mfem::BackwardEulerSolver;
  this->problem->ode_solver->Init(*(this->problem->td_operator));
}

} // namespace hephaestus

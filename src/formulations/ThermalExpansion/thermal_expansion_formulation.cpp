#include "thermal_expansion_formulation.hpp"

namespace hephaestus {

ThermalExpansionFormulation::ThermalExpansionFormulation() : SteadyStateFormulation() {

  temp_var_name = std::string("temperature");
  displacement_var_name = std::string("displacement");
  lame_param_coef_name = std::string("lame_param");
  shear_modulus_coef_name = std::string("shear_modulus");
  thermal_expansion_coef_name = std::string("thermal_expansion_coef");
  thermal_conductivity_coef_name = std::string("thermal_conductivity");
  stress_free_temp_coef_name = std::string("stress_free_temp");
  thermal_expansion_bilin_coef_name = std::string("thermal_expansion_bilin_coef");
  thermal_expansion_lin_coef_name = std::string("thermal_expansion_lin_coef");
  zero_coef_name = std::string("zero");
}

void ThermalExpansionFormulation::RegisterCoefficients() {
  std::cout << "REGISTER" << std::endl;
  hephaestus::Coefficients &coefficients = this->GetProblem()->coefficients;

  if (!coefficients.scalars.Has(lame_param_coef_name)) {
    MFEM_ABORT(lame_param_coef_name + " coefficient not found.");
  }

  if (!coefficients.scalars.Has(shear_modulus_coef_name)) {
    MFEM_ABORT(shear_modulus_coef_name + " coefficient not found.");
  }

  if (!coefficients.scalars.Has(stress_free_temp_coef_name)) {
    MFEM_ABORT(stress_free_temp_coef_name + " coefficient not found.");
  }

  if (!coefficients.scalars.Has(thermal_expansion_coef_name)) {
    MFEM_ABORT(thermal_expansion_coef_name + " coefficient not found.");
  }

  if (!coefficients.scalars.Has(thermal_conductivity_coef_name)) {
    MFEM_ABORT(thermal_conductivity_coef_name + " coefficient not found.");
  }

  if (!coefficients.scalars.Has(zero_coef_name)) {
    MFEM_ABORT(zero_coef_name + " coefficient not found.");
  }

  if (!coefficients.scalars.Has("materialTerm")) {
    mfem::SumCoefficient* materialTerm = new mfem::SumCoefficient(*(coefficients.scalars.Get(lame_param_coef_name)), *(coefficients.scalars.Get(shear_modulus_coef_name)), 3, 2);
    coefficients.scalars.Register("materialTerm", materialTerm, true);
  }

  if (!coefficients.scalars.Has("thexpStressFreeTemp")) {
    mfem::SumCoefficient* thexpStressFreeTemp = new mfem::SumCoefficient(*(coefficients.scalars.Get(thermal_expansion_coef_name)), *(coefficients.scalars.Get(stress_free_temp_coef_name)));
    coefficients.scalars.Register("thexpStressFreeTemp", thexpStressFreeTemp, true);
  }

  if (!coefficients.scalars.Has("bilinearFormCoefPositive")) {
    mfem::ProductCoefficient *bilinearFormCoefPositive = new mfem::ProductCoefficient(*(coefficients.scalars.Get(thermal_expansion_coef_name)),
                                                    *(coefficients.scalars.Get("materialTerm")));
    coefficients.scalars.Register("bilinearFormCoefPositive", bilinearFormCoefPositive, true);
  }
                                            
  if (!coefficients.scalars.Has(thermal_expansion_bilin_coef_name)) {
    coefficients.scalars.Register(thermal_expansion_bilin_coef_name,
                                    new mfem::ProductCoefficient(-1.0, *(coefficients.scalars.Get("bilinearFormCoefPositive"))), true);
  }

  if (!coefficients.scalars.Has(thermal_expansion_lin_coef_name)) {
    coefficients.scalars.Register(thermal_expansion_lin_coef_name,
                                    new mfem::ProductCoefficient(*(coefficients.scalars.Get("thexpStressFreeTemp")),
                                    *(coefficients.scalars.Get("materialTerm"))), true);
  }
}

void ThermalExpansionFormulation::ConstructEquationSystem() {
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("TempVarName", temp_var_name);
  weak_form_params.SetParam("DisplacementVarName", displacement_var_name);
  weak_form_params.SetParam("StressFreeTempCoefName", stress_free_temp_coef_name);
  weak_form_params.SetParam("LameParamCoefName", lame_param_coef_name);
  weak_form_params.SetParam("ShearModulusCoefName", shear_modulus_coef_name);
  weak_form_params.SetParam("ThermalExpansionCoefName", thermal_expansion_coef_name);
  weak_form_params.SetParam("ThermalConductivityCoefName", thermal_conductivity_coef_name);
  weak_form_params.SetParam("ThermalExpansionBilinCoefName", thermal_expansion_bilin_coef_name);
  weak_form_params.SetParam("ThermalExpansionLinCoefName", thermal_expansion_lin_coef_name);
  weak_form_params.SetParam("ZeroCoefName", zero_coef_name);
  this->GetProblem()->eq_sys =
      std::make_unique<hephaestus::ThermalExpansionEquationSystem>(weak_form_params);
}


void ThermalExpansionFormulation::ConstructOperator() {
  hephaestus::InputParameters &solver_options =

      this->GetProblem()->solver_options;
  solver_options.SetParam("TempVarName", temp_var_name);
  solver_options.SetParam("DisplacementVarName", displacement_var_name);
  solver_options.SetParam("LameCoefName", lame_param_coef_name);
  solver_options.SetParam("ShearModulusCoefName", shear_modulus_coef_name);
  solver_options.SetParam("ThermalExpansionCoefName", thermal_expansion_coef_name);
  solver_options.SetParam("ThermalConductivityCoefName", thermal_conductivity_coef_name);
  solver_options.SetParam("StressFreeTempCoefName", stress_free_temp_coef_name);
  solver_options.SetParam("ZeroCoefName", zero_coef_name);

  this->problem->eq_sys_operator = std::make_unique<hephaestus::ThermalExpansionOperator>(
      *(this->problem->pmesh), this->problem->fespaces,
      this->problem->gridfunctions, this->problem->bc_map,
      this->problem->coefficients, this->problem->sources,
      this->problem->solver_options);
  this->problem->eq_sys_operator->SetEquationSystem(
      this->problem->eq_sys.get());
  this->problem->eq_sys_operator->SetGridFunctions();
};

void ThermalExpansionFormulation::RegisterGridFunctions() {
  int &myid = this->GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = this->GetProblem()->gridfunctions;
  hephaestus::FESpaces &fespaces = this->GetProblem()->fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(temp_var_name)) {
    if (myid == 0) {
      MFEM_WARNING(temp_var_name << " not found in gridfunctions: building "
                                      "gridfunction from defaults");
    }
    AddFESpace(std::string("_TempFESpace"), std::string("H1_"));
    AddGridFunction(temp_var_name, std::string("_TempFESpace"));
  };

  if (!gridfunctions.Has(displacement_var_name)) {
      if (myid == 0) {
        MFEM_WARNING(displacement_var_name << " not found in gridfunctions: building "
                                        "gridfunction from defaults");
      }
      AddFESpace(std::string("_DisplacementFESpace"), std::string("H1_"), 3, mfem::Ordering::byVDIM);
      AddGridFunction(displacement_var_name, std::string("_DisplacementFESpace"));
  };
    // Register time derivatives
  SteadyStateProblemBuilder::RegisterGridFunctions();
};



/**
 * OPERATOR
 * 
 * 
 * 
 * 
 * 
 * 
*/
ThermalExpansionOperator::ThermalExpansionOperator(
    mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
    hephaestus::GridFunctions &gridfunctions, hephaestus::BCMap &bc_map,
    hephaestus::Coefficients &coefficients, hephaestus::Sources &sources,
    hephaestus::InputParameters &solver_options)    
    : EquationSystemOperator(pmesh, fespaces, gridfunctions,
                                            bc_map, coefficients, sources,
                                            solver_options)
      {}



// This method no longer needs to be overridden now that eq is being used
// void ThermalExpansionOperator::SetGridFunctions() {
//   state_var_names.push_back(temp_var_name);
//   state_var_names.push_back(displacement_var_name);
  
//   EquationSystemOperator::SetGridFunctions();

//   t_ = new mfem::ParGridFunction(local_test_vars.at(0)->ParFESpace());
//   u_ = new mfem::ParGridFunction(local_test_vars.at(1)->ParFESpace());
// };

void ThermalExpansionOperator::Init(mfem::Vector &X) {
  EquationSystemOperator::Init(X);
}

void ThermalExpansionOperator::Solve(mfem::Vector &X) {  

  _equation_system->FormLinearSystem(blockA, trueX, trueRhs);

  preconditioner = new mfem::HypreBoomerAMG(*blockA.As<mfem::HypreParMatrix>());
  preconditioner->SetSystemsOptions(this->pmesh_->Dimension());
  solver = new mfem::HyprePCG(*blockA.As<mfem::HypreParMatrix>());
  
  solver->SetPreconditioner(*preconditioner);
  solver->SetTol(1e-11);
  solver->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _gridfunctions);
}

/**
 * EQUATION SYSTEM
 * 
 * 
 * 
 * 
 * 
 * 
*/
ThermalExpansionEquationSystem::ThermalExpansionEquationSystem(
  const hephaestus::InputParameters &params)
  : EquationSystem(params),
    temp_var_name(params.GetParam<std::string>("TempVarName")), 
    displacement_var_name(params.GetParam<std::string>("DisplacementVarName")), 
    stress_free_temp_coef_name(params.GetParam<std::string>("StressFreeTempCoefName")),
    thermal_conductivity_coef_name(params.GetParam<std::string>("ThermalConductivityCoefName")),
    lame_param_coef_name(params.GetParam<std::string>("LameParamCoefName")),
    shear_modulus_coef_name(params.GetParam<std::string>("ShearModulusCoefName")),
    thermal_expansion_bilin_coef_name(params.GetParam<std::string>("ThermalExpansionBilinCoefName")),
    thermal_expansion_lin_coef_name(params.GetParam<std::string>("ThermalExpansionLinCoefName"))
     {}

void ThermalExpansionEquationSystem::Init(hephaestus::GridFunctions &gridfunctions,
                                  const hephaestus::FESpaces &fespaces,
                                  hephaestus::BCMap &bc_map,
                                  hephaestus::Coefficients &coefficients) {

  EquationSystem::Init(gridfunctions, fespaces, bc_map, coefficients);
}

void ThermalExpansionEquationSystem::addKernels() {
  // Add missing variable names
  addVariableNameIfMissing(temp_var_name);
  addVariableNameIfMissing(displacement_var_name);


  hephaestus::InputParameters diffusionIntegratorParams;
  diffusionIntegratorParams.SetParam("CoefficientName", thermal_conductivity_coef_name);
  addKernel(temp_var_name, new hephaestus::DiffusionKernel(diffusionIntegratorParams));

  hephaestus::InputParameters tempLFParams;
  tempLFParams.SetParam("CoefficientName", zero_coef_name);
  addKernel(temp_var_name, new hephaestus::DomainLFKernel(tempLFParams));

  hephaestus::InputParameters elasticityIntegratorParams;
  elasticityIntegratorParams.SetParam("LameParameterCoefName", lame_param_coef_name);
  elasticityIntegratorParams.SetParam("ShearModulusCoefName", shear_modulus_coef_name);
  addKernel(displacement_var_name, new hephaestus::LinearElasticityKernel(elasticityIntegratorParams));

  hephaestus::InputParameters domainDivergenceLFParams;
  domainDivergenceLFParams.SetParam("CoefficientName", thermal_expansion_lin_coef_name);
  addKernel(displacement_var_name, new hephaestus::DomainDivergenceLFKernel(domainDivergenceLFParams));

  hephaestus::InputParameters mixedWeakDivergenceParams;
  mixedWeakDivergenceParams.SetParam("CoefficientName", thermal_expansion_bilin_coef_name);
  addKernel(temp_var_name, displacement_var_name,
            new hephaestus::MixedWeakDivergenceKernel(mixedWeakDivergenceParams));
}

} //hephaestus

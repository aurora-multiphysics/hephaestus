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
}

void ThermalExpansionFormulation::RegisterCoefficients() {
  
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

  mfem::SumCoefficient materialTerm(*(coefficients.scalars.Get(lame_param_coef_name)),
                                     *(coefficients.scalars.Get(shear_modulus_coef_name)), 3, 2);

  mfem::SumCoefficient thexpStressFreeTemp(*(coefficients.scalars.Get(thermal_expansion_coef_name)),
                                            *(coefficients.scalars.Get(stress_free_temp_coef_name)));

  mfem::ProductCoefficient bilinearFormCoefPositive(*(coefficients.scalars.Get(thermal_expansion_coef_name)),
                                                    materialTerm);

  if (!coefficients.scalars.Has(thermal_expansion_bilin_coef_name)) {
    coefficients.scalars.Register(thermal_expansion_bilin_coef_name,
                                    new mfem::ProductCoefficient(-1, bilinearFormCoefPositive), true);
  }

  if (!coefficients.scalars.Has(thermal_expansion_lin_coef_name)) {
    coefficients.scalars.Register(thermal_expansion_lin_coef_name,
                                    new mfem::ProductCoefficient(*(coefficients.scalars.Get(stress_free_temp_coef_name)), materialTerm), true);
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
                                            solver_options),
      temp_var_name(solver_options.GetParam<std::string>("TempVarName")),
      displacement_var_name(solver_options.GetParam<std::string>("DisplacementVarName")),
      lame_coef_name(solver_options.GetParam<std::string>("LameCoefName")),
      shear_modulus_coef_name(solver_options.GetParam<std::string>("ShearModulusCoefName")),
      thermal_expansion_coef_name(solver_options.GetParam<std::string>("ThermalExpansionCoefName")), 
      thermal_conductivity_coef_name(solver_options.GetParam<std::string>("ThermalConductivityCoefName")),
      stress_free_temp_coef_name(solver_options.GetParam<std::string>("StressFreeTempCoefName")) {}



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
  lameCoef_ = _coefficients.scalars.Get(lame_coef_name);
  shearModulusCoef_ = _coefficients.scalars.Get(shear_modulus_coef_name);
  thermalExpansionCoef_ = _coefficients.scalars.Get(thermal_expansion_coef_name);
  stressFreeTempCoef_ = _coefficients.scalars.Get(stress_free_temp_coef_name);
  thermalConductivityCoef_ = _coefficients.scalars.Get(thermal_conductivity_coef_name);
}

void ThermalExpansionOperator::Solve(mfem::Vector &X) {  

  _equation_system->FormLinearSystem(blockA, trueX, trueRhs);

  preconditioner = new mfem::HypreBoomerAMG(*blockA.As<mfem::HypreParMatrix>());
  solver = new mfem::HyprePCG(MPI_COMM_WORLD);
  solver->SetOperator(*blockA.As<mfem::HypreParMatrix>());
  solver->SetPreconditioner(*preconditioner);
  solver->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _gridfunctions);

//   solver.SetOperator(*A1);
//   solver.SetPreconditioner(*amg);
//   solver.SetPrintLevel(2);
//   solver.Mult(trueRhs, trueX);
}

// void ThermalExpansionOperator::Solve(mfem::Vector &X) {  

//   mfem::Array2D<mfem::HypreParMatrix *> OpBlocks;
//   OpBlocks.DeleteAll();
//   OpBlocks.SetSize(2,2);
//   OpBlocks(0, 0) = new mfem::HypreParMatrix;
//   OpBlocks(1, 1) = new mfem::HypreParMatrix;
//   OpBlocks(1, 0) = new mfem::HypreParMatrix;
//   OpBlocks(0, 1) = nullptr;


//   mfem::Array<int> offsets({0, t_->ParFESpace()->TrueVSize(), 
//                            t_->ParFESpace()->TrueVSize() + u_->ParFESpace()->TrueVSize()});

//   // Apply dirichlet BC's to temperature and displacement grid functions
//   _bc_map.applyEssentialBCs(temp_var_name, ess_temp_tdofs_, *t_,
//                             pmesh_);
                            
//   _bc_map.applyEssentialBCs(displacement_var_name, ess_disp_tdofs_, *u_,
//                             pmesh_);

//   // Set up bilinear forms
//   aMixed_ = new mfem::ParMixedBilinearForm(t_->ParFESpace(), u_->ParFESpace());
//   a1_ = new mfem::ParBilinearForm(t_->ParFESpace());
//   a2_ = new mfem::ParBilinearForm(u_->ParFESpace());
//   // Set up linear forms
//   b1_ = new mfem::ParLinearForm(t_->ParFESpace());
//   b2_ = new mfem::ParLinearForm(u_->ParFESpace());


//   // Manipulate existing coefficients to get terms needed for thermal expansion.
//   // Bilinear Form Coef: - α * (3λ + 2μ)
//   // Linear Form Coef:  T_{stress free} * α * (3λ + 2μ)
//   mfem::SumCoefficient materialTerm(*lameCoef_, *shearModulusCoef_, 3, 2);
//   mfem::SumCoefficient thexpStressFreeTemp(*thermalExpansionCoef_, *stressFreeTempCoef_);
//   mfem::ProductCoefficient bilinearFormCoefPositive(*thermalExpansionCoef_, materialTerm);
//   bilinearFormCoef_ = new mfem::ProductCoefficient(-1, bilinearFormCoefPositive); 
//   linearFormCoef_ = new mfem::ProductCoefficient(*stressFreeTempCoef_, materialTerm);
  
//   a1_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*thermalConductivityCoef_));
//   a2_->AddDomainIntegrator(new mfem::ElasticityIntegrator(*lameCoef_, *shearModulusCoef_));
//   aMixed_->AddDomainIntegrator(new mfem::MixedWeakDivergenceIntegrator(*bilinearFormCoef_));
  
//   b2_->AddDomainIntegrator(new mfem::DomainLFH1DivIntegrator(*linearFormCoef_));

//   a1_->Assemble();
//   a1_->Finalize();

//   a2_->Assemble();
//   a2_->Finalize();

//   aMixed_->Assemble();
//   aMixed_->Finalize();

//   b1_->Assemble();
//   b2_->Assemble();
 
//   a1_->FormLinearSystem(ess_temp_tdofs_, *t_, *b1_, *OpBlocks(0, 0), trueX.GetBlock(0), trueRhs.GetBlock(0));
//   a2_->FormLinearSystem(ess_disp_tdofs_, *u_, *b2_, *OpBlocks(1, 1), trueX.GetBlock(1), trueRhs.GetBlock(1));
//   aMixed_->FormRectangularLinearSystem(ess_temp_tdofs_, ess_disp_tdofs_, *t_, *b2_, *OpBlocks(1, 0), trueX.GetBlock(0), trueRhs.GetBlock(1));

//   mfem::HypreParMatrix *A1 = mfem::HypreParMatrixFromBlocks(OpBlocks);
//   mfem::HypreBoomerAMG *amg = new mfem::HypreBoomerAMG(*A1);
//   mfem::HyprePCG solver(MPI_COMM_WORLD);
//   solver.SetOperator(*A1);
//   solver.SetPreconditioner(*amg);
//   solver.SetPrintLevel(2);
//   solver.Mult(trueRhs, trueX);

//   delete(amg);
//   // delete(A1);
//   OpBlocks.DeleteAll();

//   a1_->RecoverFEMSolution(trueX.GetBlock(0), *b1_, *t_);
//   a2_->RecoverFEMSolution(trueX.GetBlock(1), *b2_, *u_);
  
//   *_gridfunctions.Get(state_var_names.at(0)) = *t_;
//   *_gridfunctions.Get(state_var_names.at(1)) = *u_;
// }



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
    thermal_expansion_lin_coef_name(params.GetParam<std::string>("ThermalExpansionLinCoefName")) {}

void ThermalExpansionEquationSystem::Init(hephaestus::GridFunctions &gridfunctions,
                                  const hephaestus::FESpaces &fespaces,
                                  hephaestus::BCMap &bc_map,
                                  hephaestus::Coefficients &coefficients) {

  EquationSystem::Init(gridfunctions, fespaces, bc_map, coefficients);
}

void ThermalExpansionEquationSystem::addKernels() {
  // Add missing variable names
  addVariableNameIfMissing(temp_var_name);
  // addVariableNameIfMissing(displacement_var_name);

  // hephaestus::InputParameters mixedWeakDivergenceParams;
  // mixedWeakDivergenceParams.SetParam("CoefficientName", thermal_expansion_bilin_coef_name);
  // addKernel(temp_var_name, displacement_var_name,
  //           new hephaestus::MixedWeakDivergenceKernel(mixedWeakDivergenceParams));

  hephaestus::InputParameters diffusionIntegratorParams;
  diffusionIntegratorParams.SetParam("CoefficientName", thermal_conductivity_coef_name);
  addKernel(temp_var_name, new hephaestus::DiffusionKernel(diffusionIntegratorParams));

  // hephaestus::InputParameters domainDivergenceLFParams;
  // domainDivergenceLFParams.SetParam("CoefficientName", thermal_expansion_lin_coef_name);
  // addKernel(temp_var_name, new hephaestus::DomainDivergenceLFKernel(domainDivergenceLFParams));

  // hephaestus::InputParameters elasticityIntegratorParams;
  // elasticityIntegratorParams.SetParam("LameParameterCoefName", lame_param_coef_name);
  // elasticityIntegratorParams.SetParam("ShearModulusCoefName", shear_modulus_coef_name);
  // addKernel(displacement_var_name, new hephaestus::LinearElasticityKernel(elasticityIntegratorParams));
}

} //hephaestus

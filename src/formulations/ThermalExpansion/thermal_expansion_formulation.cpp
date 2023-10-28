#include "thermal_expansion_formulation.hpp"

namespace hephaestus {

ThermalExpansionFormulation::ThermalExpansionFormulation() : TimeDomainFormulation() {

  temp_var_name = std::string("temp_var");
  displacement_var_name = std::string("displacement_var");
  stress_free_temp_coef_name = std::string("stress_free_temp");
  lame_param_coef_name = std::string("lame_param");
  shear_modulus_coef_name = std::string("shear_modulus");
  thermal_expansion_coef_name = std::string("thermal_expansion");
  thermal_expansion_coef_name = std::string("thermal_conductivity");
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
}

// void ThermalExpansionFormulation::ConstructEquationSystem() {
//   hephaestus::InputParameters weak_form_params;
//   weak_form_params.SetParam("TempVarName", temp_var_name);
//   weak_form_params.SetParam("DisplacementVarName", displacement_var_name);
//   weak_form_params.SetParam("StressFreeTempVarName", stress_free_temp_coef_name);
//   weak_form_params.SetParam("LameParamCoefName", lame_param_coef_name);
//   weak_form_params.SetParam("ShearModulusCoefName", shear_modulus_coef_name);
//   weak_form_params.SetParam("ThermalExpansionCoefName", thermal_expansion_coef_name);
//   this->GetProblem()->ss_equation_system =
//       std::make_unique<hephaestus::CurlCurlEquationSystem>(weak_form_params);
// }

void ThermalExpansionFormulation::ConstructOperator() {
  this->problem->ss_operator = std::make_unique<hephaestus::ThermalExpansionOperator>(
      *(this->problem->pmesh), this->problem->fespaces,
      this->problem->gridfunctions, this->problem->bc_map,
      this->problem->coefficients, this->problem->sources,
      this->problem->solver_options);
  this->problem->ss_operator->SetEquationSystem(
      this->problem->ss_equation_system.get());
  this->problem->ss_operator->SetGridFunctions();
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
      AddFESpace(std::string("_DisplacementFESpace"), std::string("H1_"), 3, mfem::ordering::byVdim);
      AddGridFunction(displacement_var_name, std::string("_DisplacementFESpace"));
    };
    // Register time derivatives
  SteadyStateProblemBuilder::RegisterGridFunctions();
};



ThermalExpansionOperator::ThermalExpansionOperator(
    mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
    hephaestus::GridFunctions &gridfunctions, hephaestus::BCMap &bc_map,
    hephaestus::Coefficients &coefficients, hephaestus::Sources &sources,
    hephaestus::InputParameters &solver_options)    
    : SteadyStateEquationSystemOperator(pmesh, fespaces, gridfunctions,
                                            bc_map, coefficients, sources,
                                            solver_options),
      temp_var_name(solver_options.GetParam<std::string>("TempVarName")),
      displacement_var_name(solver_options.GetParam<std::string>("DisplacementVarName")),
      lame_coef_name(solver_options.GetParam<std::string>("LameCoefName")),
      shear_modulus_coef_name(solver_options.GetParam<std::string>("ShearModulusCoefName"))
      thermal_expansion_coef_name(solver_options.GetParam<std::string>("ThermalExpansionCoefName")) 
      stress_free_temp_coef_name(solver_options.GetParam<std::string>("StressFreeTemp")) {}


void ThermalExpansionOperator::SetGridFunctions() {
  state_var_names.push_back(temp_var_name);
  state_var_names.push_back(displacement_var_name);
  
  SteadyStateEquationSystemOperator::SetGridFunctions();

  t_ = new mfem::ParGridFunction(local_test_vars.at(0)->ParFESpace());
  u_ = new mfem::ParGridFunction(local_test_vars.at(1)->ParFESpace());
};

void ThermalExpansionOperator::Init(mfem::Vector &X) {
   ::Init(X);

  lameCoef_ = _coefficients.scalars.Get(lame_coef_name);
  shearModulusCoef_ = _coefficients.scalars.Get(shear_modulus_coef_name);
  thermalExpansionCoef_ = _coefficients.scalars.Get(thermal_expansion_coef_name);
  stressFreeTempCoef_ = _coefficients.scalars.Get(stress_free_temp_coef_name);
  thermalConductivityCoef_ = _coefficients.scalars.Get(thermal_conductivity_coef_name);

}

void ThermalExpansionOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Array2D<mfem::HypreParMatrix *> hBlocks;
  hBlocks.DeleteAll();
  hBlocks.SetSize(2,2);
  hBlocks(0, 0) = new mfem::HypreParMatrix;
  hBlocks(1, 1) = new mfem::HypreParMatrix;
  hBlocks(1, 0) = new mfem::HypreParMatrix;
  hBlocks(0, 1) = nullptr;

  mfem::Array<int> offsets({0, t_->ParFESpace().TrueVSize(), 
                           t_->ParFESpace().TrueVSize() + u_->ParFESpace().TrueVSize()});

  mfem::BlockVector trueX(offsets);
  mfem::BlockVector trueRHS(offsets);

  // Set up bilinear forms
  aMixed_ = new mfem::ParMixedBilinearForm(t_->ParFESpace(), u_->ParFESpace());
  a1_ = new mfem::ParBilinearForm(t_->ParFESpace());
  a2_ = new mfem::ParBilinearForm(u_->ParFESpace());

  if(thermalExpansionCoef_) {
    aMixed_->AddDomainIntegrator(new mfem::MixedWeakDivIntegrator(lameCoef_));
  }

  // 
  a1_->AddDomainIntegrator(new mfem::DiffusionIntegrator(thermalConductivityCoef_));
  a2_->AddDomainIntegrator(new mfem::ElasticityIntegrator(lameCoef_, shearModulusCoef_));

  _bc_map.applyEssentialBCs(std::string(temp_var_name), ess_temp_bdr_tdofs_, *t_,
                            pmesh_);

  _bc_map.applyEssentialBCs(std::string(displacement_var_name), ess_disp_bdr_tdofs_, *u_,
                            pmesh_);


  _bc_map.applyIntegratedBCs(std::string(temp_var_name), *a1_, pmesh_);
  _bc_map.applyIntegratedBCs(std::string(displacement_var_name), *a2_, pmesh_);
  

  a1_->Assemble();
  a1_->Finalize();

  a2_->Assemble();
  a2_->Finalize();

  aMixed_->Assemble();
  aMixed_->Finalize();

  // Set up linear forms
  b1_ = new mfem::ParLinearForm(t_->ParFESpace());
  b2_ = new mfem::ParLinearForm(u_->ParFESpace());

  b1_->Assemble();
  b2_->Assemble();
  
  a1_->FormLinearSystem(ess_temp_bdr_tdofs, *t_, *b1_, *hBlocks(0, 0), trueX.GetBlock(0) trueRHS.GetBlock(0));
  a2_->FormLinearSystem(ess_disp_bdr_tdofs_, *u_, *b1_, *hBlocks(1, 1), trueX.GetBlock(1), trueRHS.GetBlock(1));
  aMixed_->FormRectangularLinearSystem(ess_temp_bdr_tdofs, ess_disp_bdr_tdofs_, t_, b2_, *hBlocks(1, 0), trueX.GetBlock(0), trueB.GetBlock(1));
  A1 = mfem::HypreParMatrixFromBlocks(hBlocks);
  mfem::HypreBoomerAMG *amg = new mfem::HypreBoomerAMG(A1);
  mfem::HyprePCG solver(MPI_COMM_WORLD);
  solver.SetOperator(A1);
  solver.Mult(trueRHS, trueX);
  delete A1;

  a1_->RecoverFEMSolution(trueX.GetBlock(0), *b1_, *u_);
  a2_->RecoverFEMSolution(trueX.GetBlock(1), *b2_, *t_);
}



} //hephaestus

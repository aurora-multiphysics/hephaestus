#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "kernels.hpp"
#include "inputs.hpp"
#include "sources.hpp"
#include "integrators.hpp"

namespace hephaestus {

class ThermalExpansionFormulation : public SteadyStateFormulation {
public:
  ThermalExpansionFormulation();

  // virtual void ConstructEquationSystem() override {};

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;

protected:
  std::string temp_var_name, displacement_var_name, 
                stress_free_temp_coef_name, lame_param_coef_name, shear_modulus_coef_name, thermal_expansion_coef_name, thermal_conductivity_coef_name,
                  thermal_expansion_bilin_coef_name, thermal_expansion_lin_coef_name;
};


// Do this later
class ThermalExpansionEquationSystem : public EquationSystem {
public:
  ThermalExpansionEquationSystem(const hephaestus::InputParameters &params);

  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void addKernels() override;

  std::string temp_var_name, displacement_var_name, 
                stress_free_temp_coef_name, lame_param_coef_name, shear_modulus_coef_name,
                thermal_conductivity_coef_name, thermal_expansion_bilin_coef_name, thermal_expansion_lin_coef_name;
};

class ThermalExpansionOperator : public EquationSystemOperator {
public:
  ThermalExpansionOperator(mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
                hephaestus::GridFunctions &gridfunctions,
                hephaestus::BCMap &bc_map,
                hephaestus::Coefficients &coefficients,
                hephaestus::Sources &sources,
                hephaestus::InputParameters &solver_options);

  ~ThermalExpansionOperator(){};

  virtual void SetGridFunctions() override;
  virtual void Init(mfem::Vector &X) override;
  virtual void Solve(mfem::Vector &X) override;

  // Method for manipulating existing coefficients to get the terms needed for the thermal expansion solve
  void MakeCoefficients();

  std::string temp_var_name, displacement_var_name, 
                stress_free_temp_coef_name, lame_coef_name, shear_modulus_coef_name,
                  thermal_expansion_coef_name, thermal_conductivity_coef_name,
                    thermal_expansion_bilin_coef_name, thermal_expansion_lin_coef_name;
  mfem::ParGridFunction *u_, *t_;
  mfem::ParLinearForm *b1_, *b2_;
  mfem::ParBilinearForm *a1_, *a2_;
  mfem::ParMixedBilinearForm *aMixed_; //  
 
  mfem::Coefficient *thermalConductivityCoef_;
  mfem::Coefficient *lameCoef_;  // Lame's first parameter
  mfem::Coefficient *shearModulusCoef_;  // Shear modulus 
  mfem::Coefficient *thermalExpansionCoef_;  // Thermal expansion coefficient
  mfem::Coefficient *stressFreeTempCoef_;  // Stress free temperature
  mfem::Coefficient *bilinearFormCoef_; //
  mfem::Coefficient *linearFormCoef_; //

  mfem::Array<int> ess_temp_tdofs_, ess_disp_tdofs_;  
}; //

} // namespace hephaestus

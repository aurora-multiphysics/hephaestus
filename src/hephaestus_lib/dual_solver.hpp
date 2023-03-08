#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"

namespace hephaestus {

class DualFormulation : public TransientFormulation {
public:
  DualFormulation();

  virtual hephaestus::TimeDependentEquationSystem *
  CreateEquationSystem() override;

  virtual hephaestus::TimeDomainEquationSystemOperator *
  CreateTimeDomainOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) override;

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) override;

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels) override{};

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;

protected:
  std::string h_curl_var_name, h_div_var_name, alpha_coef_name, beta_coef_name;
};

// class DualOperator : public TimeDomainEquationSystemOperator {
// public:
//   DualOperator(mfem::ParMesh &pmesh, int order,
//                mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//                mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//                hephaestus::BCMap &bc_map,
//                hephaestus::DomainProperties &domain_properties,
//                hephaestus::Sources &sources,
//                hephaestus::InputParameters &solver_options);

//   ~DualOperator(){};

//   void ImplicitSolve(const double dt, const mfem::Vector &X,
//                      mfem::Vector &dX_dt) override;
// };

class DualOperator : public TimeDomainEquationSystemOperator {
  virtual void
  SetMaterialCoefficients(hephaestus::DomainProperties &domain_properties);
  virtual void RegisterVariables();

public:
  DualOperator(mfem::ParMesh &pmesh, int order,
               mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
               mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
               hephaestus::BCMap &bc_map,
               hephaestus::DomainProperties &domain_properties,
               hephaestus::Sources &sources,
               hephaestus::InputParameters &solver_options);

  ~DualOperator(){};

  void Init(mfem::Vector &X) override;

  void buildA1(mfem::Coefficient *sigma, mfem::Coefficient *dtMuInv);
  void buildM1(mfem::Coefficient *sigma);
  void buildCurl(mfem::Coefficient *muInv);
  void buildGrad();

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;

  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;
  mfem::ParFiniteElementSpace *HDivFESpace_;

  std::string h_curl_var_name, h_div_var_name;
  std::string u_display_name, v_display_name;

  mfem::ParGridFunction u_, du_; // HCurl vector field
  mfem::ParGridFunction v_, dv_; // HDiv vector field

  // Sockets used to communicate with GLVis
  std::map<std::string, mfem::socketstream *> socks_;

protected:
  mfem::ParBilinearForm *a1;
  mfem::HypreParMatrix *A1;
  mfem::Vector *X1, *B1;

  mfem::ParDiscreteLinearOperator *curl;
  mfem::ParMixedBilinearForm *weakCurl;

  // temporary work vectors
  mfem::ParLinearForm *b1;

  double dt_A1;
  mfem::ConstantCoefficient dtCoef;  // Coefficient for timestep scaling
  mfem::ConstantCoefficient oneCoef; // Auxiliary coefficient
  mfem::Coefficient *alphaCoef;      // Reluctivity Coefficient
  mfem::Coefficient *dtAlphaCoef;
  mfem::Coefficient *betaCoef; // Electric Conductivity Coefficient
};
} // namespace hephaestus

#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class HCurlFormulation : public TimeDomainFormulation {
public:
  HCurlFormulation();

  virtual std::unique_ptr<hephaestus::TimeDependentEquationSystem>
  CreateTimeDependentEquationSystem() const override;

  virtual std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator>
  CreateTimeDomainEquationSystemOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) const override;

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) override;

  virtual void RegisterCoefficients() override;

protected:
  std::string h_curl_var_name, alpha_coef_name, beta_coef_name;
};

class HCurlOperator : public TimeDomainEquationSystemOperator {
public:
  HCurlOperator(mfem::ParMesh &pmesh,
                mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
                mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                hephaestus::BCMap &bc_map,
                hephaestus::DomainProperties &domain_properties,
                hephaestus::Sources &sources,
                hephaestus::InputParameters &solver_options);

  ~HCurlOperator(){};

  void ImplicitSolve(const double dt, const mfem::Vector &X,
                     mfem::Vector &dX_dt) override;
};

// class HCurlSolver : public TimeDomainFormulation {
//   virtual void SetMaterialCoefficients(
//       hephaestus::DomainProperties &domain_properties) override;
//   virtual void SetEquationSystem() override;

// public:
//   HCurlSolver(mfem::ParMesh &pmesh,
//               mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//               mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//               hephaestus::BCMap &bc_map,
//               hephaestus::DomainProperties &domain_properties,
//               hephaestus::Sources &sources,
//               hephaestus::InputParameters &solver_options);

//   ~HCurlSolver(){};

//   virtual void RegisterMissingVariables() override;
//   void ImplicitSolve(const double dt, const mfem::Vector &X,
//                      mfem::Vector &dX_dt) override;

//   std::string alpha_coef_name, beta_coef_name;

// protected:
// };
} // namespace hephaestus

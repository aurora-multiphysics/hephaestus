#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

class HCurlFormulation : public TimeDomainFormulation {
public:
  HCurlFormulation();

  virtual void ConstructEquationSystem() override;

  virtual std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator>
  CreateTimeDomainEquationSystemOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) const override;

  virtual void RegisterGridFunctions() override;

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

} // namespace hephaestus

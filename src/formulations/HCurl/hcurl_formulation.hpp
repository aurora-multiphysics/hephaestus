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

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;

protected:
  std::string h_curl_var_name, alpha_coef_name, beta_coef_name;
};

class CurlCurlEquationSystem : public TimeDependentEquationSystem {
public:
  CurlCurlEquationSystem(const hephaestus::InputParameters &params);

  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void addKernels() override;

  std::string h_curl_var_name, alpha_coef_name, beta_coef_name,
      dtalpha_coef_name;
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

#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "materials.hpp"

namespace hephaestus {

template <typename T> class Kernel {
public:
  Kernel() {}
  Kernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties){};

  std::string variable_name;
  virtual void Apply(T *form) = 0;
};

/*
(α∇×u_{n}, ∇×u')
*/
class WeakCurlCurlKernel : public Kernel<mfem::ParLinearForm> {
public:
  WeakCurlCurlKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParLinearForm *lf) override;

  std::string gf_name, coef_name;
  mfem::ParGridFunction *u_; //
  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::Coefficient *coef;
  mfem::ParBilinearForm *curlCurl;
};

/*
(α∇×u, ∇×u')
*/
class CurlCurlKernel : public Kernel<mfem::ParBilinearForm> {
public:
  CurlCurlKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

/*
(βu, u')
*/
class VectorFEMassKernel : public Kernel<mfem::ParBilinearForm> {
public:
  VectorFEMassKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

/*
(σ ∇ V, ∇ V')
*/
class DiffusionKernel : public Kernel<mfem::ParBilinearForm> {
public:
  DiffusionKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

/*
(σ ∇ V, u')
*/
class MixedVectorGradientKernel : public Kernel<mfem::ParMixedBilinearForm> {
public:
  MixedVectorGradientKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParMixedBilinearForm *mblf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

/*
(σ u, ∇ V')
*/
class VectorFEWeakDivergenceKernel : public Kernel<mfem::ParMixedBilinearForm> {
public:
  VectorFEWeakDivergenceKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParMixedBilinearForm *mblf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus

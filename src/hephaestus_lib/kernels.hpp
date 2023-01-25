#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "materials.hpp"

namespace hephaestus {

class Kernel {
public:
  Kernel() {}
  Kernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties){};

  std::string variable_name;
  virtual void ApplyKernel(mfem::ParLinearForm *lf){};
  virtual void ApplyKernel(mfem::ParBilinearForm *blf){};
};

class Kernels : public mfem::NamedFieldsMap<hephaestus::Kernel> {
public:
  void Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
            hephaestus::BCMap &bc_map,
            hephaestus::DomainProperties &domain_properties) {
    for (const auto &[name, kernel] : GetMap()) {
      kernel->Init(variables, fespaces, bc_map, domain_properties);
    }
  };

  void ApplyKernels(const std::string &variable_name_,
                    mfem::ParBilinearForm *blf) {
    for (const auto &[name, kernel] : GetMap()) {
      if (kernel->variable_name == variable_name_) {
        kernel->ApplyKernel(blf);
      }
    }
  };
  void ApplyKernels(const std::string &variable_name_,
                    mfem::ParLinearForm *lf) {
    for (const auto &[name, kernel] : GetMap()) {
      if (kernel->variable_name == variable_name_) {
        kernel->ApplyKernel(lf);
      }
    }
  };
  void ApplyKernels(mfem::ParLinearForm *lf) {
    for (const auto &[name, kernel] : GetMap()) {
      kernel->ApplyKernel(lf);
    }
  };
  void ApplyKernels(mfem::ParBilinearForm *blf) {
    for (const auto &[name, kernel] : GetMap()) {
      kernel->ApplyKernel(blf);
    }
  };
};

// class VectorFEMassKernel : public Kernel {
// public:
//   VectorFEMassKernel(const hephaestus::InputParameters &params);
//   mfem::ParBilinearForm *curlCurl;
//   virtual void
//   Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
//        const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
//        hephaestus::BCMap &bc_map,
//        hephaestus::DomainProperties &domain_properties) override {

//     curlCurl = new mfem::ParBilinearForm(HCurlFESpace_);
//     curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*MuInv));
//     curlCurl->Assemble();
//   };

//   virtual void ApplyKernel(mfem::ParBilinearForm *blf) {
//     blf->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
//   };
// };

class WeakCurlCurlKernel : public Kernel {
  // (α∇×u_{n}, ∇×u')
public:
  WeakCurlCurlKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void ApplyKernel(mfem::ParLinearForm *lf) override;

  std::string gf_name, coef_name;
  mfem::ParGridFunction *u_; //
  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::Coefficient *coef;
  mfem::ParBilinearForm *curlCurl;
};

}; // namespace hephaestus

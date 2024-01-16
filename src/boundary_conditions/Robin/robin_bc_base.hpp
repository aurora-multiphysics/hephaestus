#pragma once
#include "integrated_bc_base.hpp"

namespace hephaestus
{

class RobinBC : public IntegratedBC
{
public:
  RobinBC(const std::string & name_,
          mfem::Array<int> bdr_attributes_,
          mfem::BilinearFormIntegrator * blfi_re_,
          mfem::LinearFormIntegrator * lfi_re_,
          mfem::BilinearFormIntegrator * blfi_im_ = nullptr,
          mfem::LinearFormIntegrator * lfi_im_ = nullptr);

  // NB: assume ownership of pointers.
  std::unique_ptr<mfem::BilinearFormIntegrator> blfi_re{nullptr};
  std::unique_ptr<mfem::BilinearFormIntegrator> blfi_im{nullptr};

  virtual void applyBC(mfem::ParBilinearForm & a);
  virtual void applyBC(mfem::ParSesquilinearForm & a);
};

} // namespace hephaestus

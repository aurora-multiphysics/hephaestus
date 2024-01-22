#pragma once
#include "boundary_condition_base.hpp"

namespace hephaestus
{

class IntegratedBC : public BoundaryCondition
{
public:
  IntegratedBC(const std::string & name_, mfem::Array<int> bdr_attributes_);
  IntegratedBC(const std::string & name_,
               mfem::Array<int> bdr_attributes_,
               mfem::LinearFormIntegrator * lfi_re_,
               mfem::LinearFormIntegrator * lfi_im_ = nullptr);

  // NB: assume ownership of pointers.
  std::unique_ptr<mfem::LinearFormIntegrator> lfi_re;
  std::unique_ptr<mfem::LinearFormIntegrator> lfi_im;

  void ApplyBC(mfem::LinearForm & b) override;
  void ApplyBC(mfem::ComplexLinearForm & b) override;
  void ApplyBC(mfem::ParComplexLinearForm & b) override;
};

} // namespace hephaestus

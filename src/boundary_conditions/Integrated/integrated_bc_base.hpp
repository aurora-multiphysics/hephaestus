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

  mfem::LinearFormIntegrator * lfi_re;
  mfem::LinearFormIntegrator * lfi_im;

  virtual void applyBC(mfem::LinearForm & b) override;
  virtual void applyBC(mfem::ComplexLinearForm & b) override;
  virtual void applyBC(mfem::ParComplexLinearForm & b) override;
};

} // namespace hephaestus

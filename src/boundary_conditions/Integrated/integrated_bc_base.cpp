#include "integrated_bc_base.hpp"

namespace hephaestus {

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_,
                           mfem::LinearFormIntegrator *lfi_re_,
                           mfem::LinearFormIntegrator *lfi_im_)
    : BoundaryCondition(name_, bdr_attributes_), lfi_re(lfi_re_),
      lfi_im(lfi_im_) {}

void IntegratedBC::applyBC(mfem::LinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re, markers);
}

void IntegratedBC::applyBC(mfem::ComplexLinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re, lfi_im, markers);
}

void IntegratedBC::applyBC(mfem::ParComplexLinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re, lfi_im, markers);
}

} // namespace hephaestus

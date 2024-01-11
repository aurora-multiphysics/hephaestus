#include "integrated_bc_base.hpp"

namespace hephaestus {

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_,
                           mfem::LinearFormIntegrator *lfi_re_,
                           mfem::LinearFormIntegrator *lfi_im_)
    : BoundaryCondition(name_, bdr_attributes_),
      lfi_re(std::unique_ptr<mfem::LinearFormIntegrator>(lfi_re_)),
      lfi_im(std::unique_ptr<mfem::LinearFormIntegrator>(lfi_im_)) {}

void IntegratedBC::applyBC(mfem::LinearForm &b) {
  // NB: release ownership to prevent double-free. LinearForm assumes ownership.
  b.AddBoundaryIntegrator(lfi_re.release(), markers);
}

void IntegratedBC::applyBC(mfem::ComplexLinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re.release(), lfi_im.release(), markers);
}

void IntegratedBC::applyBC(mfem::ParComplexLinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re.release(), lfi_im.release(), markers);
}

} // namespace hephaestus

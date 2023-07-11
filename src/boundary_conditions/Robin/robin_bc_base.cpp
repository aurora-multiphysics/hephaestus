#include "robin_bc_base.hpp"

namespace hephaestus {

RobinBC::RobinBC(const std::string &name_, mfem::Array<int> bdr_attributes_,
                 mfem::BilinearFormIntegrator *blfi_re_,
                 mfem::LinearFormIntegrator *lfi_re_,
                 mfem::BilinearFormIntegrator *blfi_im_,
                 mfem::LinearFormIntegrator *lfi_im_)
    : IntegratedBC(name_, bdr_attributes_, lfi_re_, lfi_im_), blfi_re(blfi_re_),
      blfi_im(blfi_im_) {}

void RobinBC::applyBC(mfem::ParBilinearForm &a) {
  a.AddBoundaryIntegrator(blfi_re, markers);
}

void RobinBC::applyBC(mfem::ParSesquilinearForm &a) {
  a.AddBoundaryIntegrator(blfi_re, blfi_im, markers);
}

} // namespace hephaestus

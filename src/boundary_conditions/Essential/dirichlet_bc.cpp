#include "dirichlet_bc.hpp"

namespace hephaestus {

DirichletBC::DirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_)
    : EssentialBC(name_, bdr_attributes_) {}

DirichletBC::DirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_,
                                         mfem::ConstantCoefficient *coeff_,
                                         mfem::ConstantCoefficient *coeff_im_)
    : EssentialBC(name_, bdr_attributes_), coeff(coeff_), coeff_im(coeff_im_) {}

void DirichletBC::applyBC(mfem::GridFunction &gridfunc,
                                  mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  gridfunc.ProjectBdrCoefficient(*(this->coeff), ess_bdrs);
}

} // namespace hephaestus

#include "function_dirichlet_bc.hpp"

namespace hephaestus {

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_)
    : EssentialBC(name_, bdr_attributes_) {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_,
                                         mfem::FunctionCoefficient *coeff_,
                                         mfem::FunctionCoefficient *coeff_im_)
    : EssentialBC(name_, bdr_attributes_), coeff(coeff_), coeff_im(coeff_im_) {}

void FunctionDirichletBC::applyBC(mfem::GridFunction &gridfunc,
                                  mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  gridfunc.ProjectBdrCoefficient(*(this->coeff), ess_bdrs);
}

} // namespace hephaestus

#include "vector_function_dirichlet_bc.hpp"

namespace hephaestus {

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_)
    : EssentialBC(name_, bdr_attributes_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_,
    mfem::VectorFunctionCoefficient *vec_coeff_,
    mfem::VectorFunctionCoefficient *vec_coeff_im_)
    : EssentialBC(name_, bdr_attributes_), vec_coeff(vec_coeff_),
      vec_coeff_im(vec_coeff_im_) {}

void VectorFunctionDirichletBC::applyBC(mfem::GridFunction &gridfunc,
                                        mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  if (this->vec_coeff == NULL) {
    MFEM_ABORT(
        "Boundary condition does not store valid coefficients to specify the "
        "components of the vector at the Dirichlet boundary.");
  }
  gridfunc.ProjectBdrCoefficientTangent(*(this->vec_coeff), ess_bdrs);
}

void VectorFunctionDirichletBC::applyBC(mfem::ParComplexGridFunction &gridfunc,
                                        mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  if (this->vec_coeff == NULL || this->vec_coeff_im == NULL) {
    MFEM_ABORT(
        "Boundary condition does not store valid coefficients to specify both "
        "the real and imaginary components of the vector at the Dirichlet "
        "boundary.");
  }
  gridfunc.ProjectBdrCoefficientTangent(*(this->vec_coeff),
                                        *(this->vec_coeff_im), ess_bdrs);
}

} // namespace hephaestus

#include "vector_dirichlet_bc.hpp"

namespace hephaestus {

VectorDirichletBC::VectorDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_)
    : EssentialBC(name_, bdr_attributes_) {}

VectorDirichletBC::VectorDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_,
    mfem::VectorCoefficient *vec_coeff_,
    mfem::VectorCoefficient *vec_coeff_im_,
    APPLY_TYPE boundary_apply_type_)
    : EssentialBC(name_, bdr_attributes_), vec_coeff(vec_coeff_),
      vec_coeff_im(vec_coeff_im_), boundary_apply_type(boundary_apply_type_) {}

void VectorDirichletBC::applyBC(mfem::GridFunction &gridfunc,
                                        mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  if (this->vec_coeff == NULL) {
    MFEM_ABORT(
        "Boundary condition does not store valid coefficients to specify the "
        "components of the vector at the Dirichlet boundary.");
  }

  if(boundary_apply_type == STANDARD)
  {
    gridfunc.ProjectBdrCoefficient(*(this->vec_coeff), ess_bdrs);
  }

  if(boundary_apply_type == TANGENTIAL)
  {
    gridfunc.ProjectBdrCoefficientTangent(*(this->vec_coeff), ess_bdrs);
  }

  if(boundary_apply_type == NORMAL)
  {
    gridfunc.ProjectBdrCoefficientNormal(*(this->vec_coeff), ess_bdrs);
  }
}

void VectorDirichletBC::applyBC(mfem::ParComplexGridFunction &gridfunc,
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

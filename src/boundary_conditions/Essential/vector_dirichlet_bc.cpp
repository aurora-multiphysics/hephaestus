#include "vector_dirichlet_bc.hpp"

namespace hephaestus
{

VectorDirichletBC::VectorDirichletBC(const std::string & name_, mfem::Array<int> bdr_attributes_)
  : EssentialBC(name_, bdr_attributes_)
{
}

VectorDirichletBC::VectorDirichletBC(const std::string & name_,
                                     mfem::Array<int> bdr_attributes_,
                                     mfem::VectorCoefficient * vec_coeff_,
                                     mfem::VectorCoefficient * vec_coeff_im_,
                                     APPLY_TYPE boundary_apply_type_)
  : EssentialBC(name_, bdr_attributes_),
    vec_coeff(vec_coeff_),
    vec_coeff_im(vec_coeff_im_),
    boundary_apply_type(boundary_apply_type_)
{
}

void
VectorDirichletBC::applyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = getMarkers(*mesh_);
  if (vec_coeff == NULL)
  {
    MFEM_ABORT("Boundary condition does not store valid coefficients to specify the "
               "components of the vector at the Dirichlet boundary.");
  }

  switch (boundary_apply_type)
  {
    case STANDARD:
      gridfunc.ProjectBdrCoefficient(*(vec_coeff), ess_bdrs);
      break;
    case NORMAL:
      gridfunc.ProjectBdrCoefficientNormal(*(vec_coeff), ess_bdrs);
      break;
    case TANGENTIAL:
      gridfunc.ProjectBdrCoefficientTangent(*(vec_coeff), ess_bdrs);
  }
}

void
VectorDirichletBC::applyBC(mfem::ParComplexGridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = getMarkers(*mesh_);
  if (vec_coeff == NULL || vec_coeff_im == NULL)
  {
    MFEM_ABORT("Boundary condition does not store valid coefficients to specify both "
               "the real and imaginary components of the vector at the Dirichlet "
               "boundary.");
  }
  gridfunc.ProjectBdrCoefficientTangent(*(vec_coeff), *(vec_coeff_im), ess_bdrs);
}

} // namespace hephaestus

#include "boundary_conditions.hpp"

namespace hephaestus {

BoundaryCondition::BoundaryCondition() {}

BoundaryCondition::BoundaryCondition(const std::string &name_,
                                     mfem::Array<int> bdr_attributes_)
    : name(name_), bdr_attributes(bdr_attributes_) {}

mfem::Array<int> BoundaryCondition::getMarkers(mfem::Mesh &mesh) {
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes,
                             markers);
  return markers;
}

FunctionDirichletBC::FunctionDirichletBC() {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_,
                                         mfem::FunctionCoefficient *coeff_,
                                         mfem::FunctionCoefficient *coeff_im_)
    : BoundaryCondition(name_, bdr_attributes_), coeff(coeff_),
      coeff_im(coeff_im_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC() {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_,
    mfem::VectorFunctionCoefficient *vec_coeff_,
    mfem::VectorFunctionCoefficient *vec_coeff_im_)
    : BoundaryCondition(name_, bdr_attributes_), vec_coeff(vec_coeff_),
      vec_coeff_im(vec_coeff_im_) {}

IntegratedBC::IntegratedBC() {}

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

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

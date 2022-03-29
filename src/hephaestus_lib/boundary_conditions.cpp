#include "boundary_conditions.hpp"

namespace hephaestus {

BoundaryCondition::BoundaryCondition() {}

BoundaryCondition::BoundaryCondition(const std::string &boundary_name,
                                     mfem::Array<int> boundary_ids)
    : name(boundary_name), bdr_attributes(boundary_ids) {}

mfem::Array<int> BoundaryCondition::getMarkers(mfem::Mesh &mesh) {
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes,
                             markers);
  return markers;
}

FunctionDirichletBC::FunctionDirichletBC() {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &boundary_name,
                                         mfem::Array<int> boundary_ids)
    : BoundaryCondition(boundary_name, boundary_ids) {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &boundary_name,
                                         mfem::Array<int> boundary_ids,
                                         mfem::FunctionCoefficient *coeff_,
                                         mfem::FunctionCoefficient *coeff_im_)
    : BoundaryCondition(boundary_name, boundary_ids), coeff(coeff_),
      coeff_im(coeff_im_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC() {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &boundary_name, mfem::Array<int> boundary_ids)
    : BoundaryCondition(boundary_name, boundary_ids) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &boundary_name, mfem::Array<int> boundary_ids,
    mfem::VectorFunctionCoefficient *vec_coeff_,
    mfem::VectorFunctionCoefficient *vec_coeff_im_)
    : BoundaryCondition(boundary_name, boundary_ids), vec_coeff(vec_coeff_),
      vec_coeff_im(vec_coeff_im_) {}

IntegratedBC::IntegratedBC() {}

IntegratedBC::IntegratedBC(const std::string &boundary_name,
                           mfem::Array<int> boundary_ids)
    : BoundaryCondition(boundary_name, boundary_ids) {}

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
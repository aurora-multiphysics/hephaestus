#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

class BoundaryCondition {
public:
  BoundaryCondition();
  BoundaryCondition(const std::string &name_, mfem::Array<int> bdr_attributes_);
  mfem::Array<int> getMarkers(mfem::Mesh &mesh);

  std::string name;
  mfem::Array<int> bdr_attributes;
  mfem::Array<int> markers;

  virtual void applyBC(mfem::LinearForm &b){};
  virtual void applyBC(mfem::ComplexLinearForm &b){};
  virtual void applyBC(mfem::ParComplexLinearForm &b){};
};

class EssentialBC : public BoundaryCondition {
public:
  EssentialBC();
  EssentialBC(const std::string &name_, mfem::Array<int> bdr_attributes_);

  virtual void applyBC(mfem::GridFunction &gridfunc, mfem::Mesh *mesh_,
                       double time = 0.0){};
};

class FunctionDirichletBC : public EssentialBC {
public:
  FunctionDirichletBC();
  FunctionDirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_);
  FunctionDirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_,
                      mfem::FunctionCoefficient *coeff_,
                      mfem::FunctionCoefficient *coeff_im_ = nullptr);

  virtual void applyBC(mfem::GridFunction &gridfunc, mfem::Mesh *mesh_,
                       double time = 0.0) override;

  mfem::FunctionCoefficient *coeff;
  mfem::FunctionCoefficient *coeff_im;
};

class VectorFunctionDirichletBC : public EssentialBC {
public:
  VectorFunctionDirichletBC();
  VectorFunctionDirichletBC(const std::string &name_,
                            mfem::Array<int> bdr_attributes_);
  VectorFunctionDirichletBC(
      const std::string &name_, mfem::Array<int> bdr_attributes_,
      mfem::VectorFunctionCoefficient *vec_coeff_,
      mfem::VectorFunctionCoefficient *vec_coeff_im_ = nullptr);

  virtual void applyBC(mfem::GridFunction &gridfunc, mfem::Mesh *mesh_,
                       double time = 0.0) override;

  mfem::VectorFunctionCoefficient *vec_coeff;
  mfem::VectorFunctionCoefficient *vec_coeff_im;
};

class IntegratedBC : public BoundaryCondition {
public:
  IntegratedBC();
  IntegratedBC(const std::string &name_, mfem::Array<int> bdr_attributes_);
  IntegratedBC(const std::string &name_, mfem::Array<int> bdr_attributes_,
               mfem::LinearFormIntegrator *lfi_re_,
               mfem::LinearFormIntegrator *lfi_im_ = nullptr);

  mfem::LinearFormIntegrator *lfi_re;
  mfem::LinearFormIntegrator *lfi_im;

  virtual void applyBC(mfem::LinearForm &b) override;
  virtual void applyBC(mfem::ComplexLinearForm &b) override;
  virtual void applyBC(mfem::ParComplexLinearForm &b) override;
};

class BCMap : public std::map<std::string, hephaestus::BoundaryCondition *> {
public:
  mfem::Array<int> getEssentialBdrMarkers(const std::string &name_,
                                          mfem::Mesh *mesh_);

  mfem::Array<int> applyEssentialBCs(const std::string &name_,
                                     mfem::GridFunction &gridfunc,
                                     mfem::Mesh *mesh_, double time = 0.0);

  void applyIntegratedBCs(const std::string &name_, mfem::LinearForm &lf,
                          mfem::Mesh *mesh_);
  void applyIntegratedBCs(const std::string &name_,
                          mfem::ParComplexLinearForm &clf, mfem::Mesh *mesh_);
};

} // namespace hephaestus

#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

class BoundaryCondition {
public:
  BoundaryCondition();
  BoundaryCondition(const std::string &boundary_name,
                    mfem::Array<int> boundary_ids);
  mfem::Array<int> getMarkers(mfem::Mesh &mesh);

  std::string name;
  mfem::Array<int> bdr_attributes;
  mfem::Array<int> markers;

  virtual void applyBC(mfem::LinearForm &b){};
  virtual void applyBC(mfem::ComplexLinearForm &b){};
  virtual void applyBC(mfem::ParComplexLinearForm &b){};
};

class FunctionDirichletBC : public BoundaryCondition {
public:
  FunctionDirichletBC();
  FunctionDirichletBC(const std::string &boundary_name,
                      mfem::Array<int> boundary_ids);
  FunctionDirichletBC(const std::string &boundary_name,
                      mfem::Array<int> boundary_ids,
                      mfem::FunctionCoefficient *coeff_,
                      mfem::FunctionCoefficient *coeff_im_ = nullptr);

  std::function<double(const mfem::Vector &, double)> scalar_func;
  std::function<double(const mfem::Vector &, double)> scalar_func_im;

  mfem::FunctionCoefficient *coeff;
  mfem::FunctionCoefficient *coeff_im;
};

class VectorFunctionDirichletBC : public BoundaryCondition {
public:
  VectorFunctionDirichletBC();
  VectorFunctionDirichletBC(const std::string &boundary_name,
                            mfem::Array<int> boundary_ids);
  VectorFunctionDirichletBC(
      const std::string &boundary_name, mfem::Array<int> boundary_ids,
      mfem::VectorFunctionCoefficient *vec_coeff_,
      mfem::VectorFunctionCoefficient *vec_coeff_im_ = nullptr);

  std::function<void(const mfem::Vector &, mfem::Vector &)> vector_func;
  std::function<void(const mfem::Vector &, mfem::Vector &)> vector_func_im;

  mfem::VectorFunctionCoefficient *vec_coeff;
  mfem::VectorFunctionCoefficient *vec_coeff_im;
};

class IntegratedBC : public BoundaryCondition {
public:
  IntegratedBC();
  IntegratedBC(const std::string &boundary_name, mfem::Array<int> boundary_ids);

  mfem::LinearFormIntegrator *lfi_re;
  mfem::LinearFormIntegrator *lfi_im;

  virtual void applyBC(mfem::LinearForm &b) override;
  virtual void applyBC(mfem::ComplexLinearForm &b) override;
  virtual void applyBC(mfem::ParComplexLinearForm &b) override;
};

class BCMap : public std::map<std::string, hephaestus::BoundaryCondition *> {
  // public:
  //   BCMap();

  //   void setBC(std::string bc_name, BoundaryCondition bc);
  //   BoundaryCondition getBC(std::string bc_name);

  //   std::map<std::string, hephaestus::BoundaryCondition *> bc_map;
};

} // namespace hephaestus
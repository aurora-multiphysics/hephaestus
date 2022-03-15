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
  std::function<double(const mfem::Vector &, double)> scalar_func;
  std::function<double(const mfem::Vector &, double)> scalar_func_im;
  std::function<void(const mfem::Vector &, mfem::Vector &)> vector_func;
  std::function<void(const mfem::Vector &, mfem::Vector &)> vector_func_im;

  mfem::Array<int> bdr_attributes;
  mfem::Array<int> markers;

  virtual void applyBC(mfem::LinearForm &b){};
  virtual void applyBC(mfem::ComplexLinearForm &b){};
  virtual void applyBC(mfem::ParComplexLinearForm &b){};
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

class BCMap {
 public:
  BCMap();

  void setBC(std::string bc_name, BoundaryCondition bc);
  BoundaryCondition getBC(std::string bc_name);

  std::map<std::string, hephaestus::BoundaryCondition *> bc_map;
};

}  // namespace hephaestus
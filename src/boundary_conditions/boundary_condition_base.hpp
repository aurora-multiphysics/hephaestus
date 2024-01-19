#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus
{

class BoundaryCondition
{
public:
  BoundaryCondition(std::string name_, mfem::Array<int> bdr_attributes_);
  mfem::Array<int> getMarkers(mfem::Mesh & mesh);

  std::string name;
  mfem::Array<int> bdr_attributes;
  mfem::Array<int> markers;

  virtual void applyBC(mfem::LinearForm & b){};
  virtual void applyBC(mfem::ComplexLinearForm & b){};
  virtual void applyBC(mfem::ParComplexLinearForm & b){};
};

} // namespace hephaestus

#pragma once
#include "mfem.hpp"
#include "named_fields_map.hpp"

namespace hephaestus
{

class FECollections : public hephaestus::NamedFieldsMap<mfem::FiniteElementCollection>
{
};

class FESpaces : public hephaestus::NamedFieldsMap<mfem::ParFiniteElementSpace>
{
public:
  ~FESpaces() override = default;

  /// @brief Update stored fespaces on mesh change.
  virtual void Update();
};

class GridFunctions : public hephaestus::NamedFieldsMap<mfem::ParGridFunction>
{
public:
  ~GridFunctions() override = default;

  /// @brief Update stored gridfunctions on mesh change.
  virtual void Update();
};

} // namespace hephaestus

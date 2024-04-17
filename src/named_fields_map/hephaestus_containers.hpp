#pragma once
#include "mfem.hpp"
#include "named_fields_map.hpp"
#include "update_interface.hpp"

namespace hephaestus
{

class FECollections : public hephaestus::NamedFieldsMap<mfem::FiniteElementCollection>
{
};

class FESpaces : public hephaestus::NamedFieldsMap<mfem::ParFiniteElementSpace>,
                 public MeshUpdateInterface
{
public:
  ~FESpaces() override = default;

  /// @brief Update stored fespaces on mesh change.
  void Update() override;
};

class GridFunctions : public hephaestus::NamedFieldsMap<mfem::ParGridFunction>,
                      public MeshUpdateInterface
{
public:
  ~GridFunctions() override = default;

  /// @brief Update stored gridfunctions on mesh change.
  void Update() override;
};

} // namespace hephaestus

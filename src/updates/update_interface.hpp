#pragma once

namespace hephaestus
{

/// @brief Simple interface with "Update" method. All classes which need to be
/// updated on a mesh change should inherit from this interface.
struct UpdateInterface
{
  virtual ~UpdateInterface() = default;

  /// @brief Implement in derived classes on mesh change.
  virtual void Update() = 0;
};

}
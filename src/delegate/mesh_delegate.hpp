#pragma once
#include <set>
#include <memory>
#include "mfem"

namespace hephaestus
{

/**
 * The delegate for mesh updates. Classes should inherit from this, implement OnMeshDidUpdate() and
 * instances should be added by calling AddDelegate().
 */
class MeshDelegate
{
public:
  MeshDelegate() = default;
  ~MeshDelegate() = default;

  virtual void OnMeshDidUpdate() { MFEM_ABORT("Not implemented."); }
};

/**
 * The delegator for mesh updates. Classes should inherit from this.
 */
class MeshDelegator
{
public:
  MeshDelegator() = default;
  ~MeshDelegator() = default;

  inline void AddDelegate(MeshDelegate * delegate) { _delegates.insert(delegate); }
  inline void RemoveDelegate(MeshDelegate * delegate) { _delegates.erase(delegate); }
  inline void ClearDelegates() { _delegates.clear(); }

  void MeshDidUpdate()
  {
    for (const auto & delegate : _delegates)
    {
      delegate->OnMeshDidUpdate();
    }
  }

private:
  std::set<MeshDelegate *> _delegates;
};
}

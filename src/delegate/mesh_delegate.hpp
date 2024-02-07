#pragma once
#include <set>
#include <memory>
#include <mfem.hpp>

namespace hephaestus
{
// Forwards declaration.
class MeshDelegate;

/**
 * The delegator for mesh updates. Classes should inherit from this.
 */
class MeshDelegator
{
public:
  MeshDelegator() = default;
  virtual ~MeshDelegator() = default;

  /// Add a new delegate. This will be called by delegates themselves. No danger if this is
  /// called multiple times (set prevents double-registration).
  inline void AddDelegate(MeshDelegate * delegate) { _delegates.insert(delegate); }

  /// Remove a delegate. This should be called in the destructor of delegates to prevent
  /// MeshDidUpdate() attempting to update a deleted delegate.
  inline void RemoveDelegate(MeshDelegate * delegate) { _delegates.erase(delegate); }

protected:
  /// NB: should be a protected method so only inherited classes can call this. We don't want
  /// delegates screwing around with their associated delegator.
  inline void ClearDelegates() { _delegates.clear(); }

  /// Called to update all delegates on a mesh update. NB: this method must be protected
  /// so that inherited classes can call this but delegates which hold a weak reference
  /// to their delegator cannot call it.
  void MeshDidUpdate();

private:
  std::set<MeshDelegate *> _delegates;
};

/**
 * The delegate for mesh updates. Classes should inherit from this, implement OnMeshDidUpdate() and
 * instances should be added by calling AddDelegate().
 */
class MeshDelegate
{
public:
  /// Delete default constructor.
  MeshDelegate() = delete;

  /// Default initializer will now ensure that we register ourselves to a mesh delegator.
  MeshDelegate(MeshDelegator & delegator) : _delegator{delegator} { _delegator.AddDelegate(this); }

  /// Destructor ensures that we are removed from the delegator.
  virtual ~MeshDelegate() { _delegator.RemoveDelegate(this); }

  virtual void OnMeshDidUpdate() { MFEM_ABORT("Not implemented."); }

private:
  /// Hold a pointer to our delegator. For now, assume that the delegator will have
  /// a lifetime at least equal to that of the longest-living delegate. We may need to
  /// implement a shared pointer to ensure that the member variable is always valid.
  MeshDelegator & _delegator;
};

/**
 * Delegate template class for wrapping around existing MFEM/Hephaestus classes. Inherits from
 * MeshDelegate.
 */
template <class T>
class MeshDelegateTemplate : public MeshDelegate, public T
{
  /// Constructor for template class.
  template <class... TArgs>
  MeshDelegateTemplate(MeshDelegator & delegator, TArgs &&... args)
    : MeshDelegate(delegator), T(std::forward<TArgs>(args)...)
  {
  }
};

// Templates.
using ParBilinearFormDelegate = MeshDelegateTemplate<mfem::ParBilinearForm>;
using ParLinearFormDelegate = MeshDelegateTemplate<mfem::ParLinearForm>;
using ParNonlinearFormDelegate = MeshDelegateTemplate<mfem::ParNonlinearForm>;
using ParMixedBilinearFormDelegate = MeshDelegateTemplate<mfem::ParMixedBilinearForm>;
}

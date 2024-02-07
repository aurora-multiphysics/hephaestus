#pragma once
#include <set>
#include <memory>
#include <mfem.hpp>

namespace hephaestus
{
// Forwards declaration.
class Listener;

/**
 * The broadcaster for mesh updates. Classes should inherit from this.
 */
class Broadcaster
{
public:
  Broadcaster() = default;
  virtual ~Broadcaster() = default;

  /// Add a new listener. This will be called by listeners themselves. No danger if this is
  /// called multiple times (set prevents double-registration).
  inline void AddListener(Listener * listener) { _listeners.insert(listener); }

  /// Remove a listener. This should be called in the destructor of listeners to prevent
  /// BroadCastUpdate() attempting to update a deleted listener.
  inline void RemoveListener(Listener * listener) { _listeners.erase(listener); }

protected:
  /// NB: should be a protected method so only inherited classes can call this. We don't want
  /// listeners screwing around with their associated broadcaster.
  inline void RemoveListeners() { _listeners.clear(); }

  /// Called to update all listeners. NB: this method must be protected so that inherited classes
  /// can call this but listeners which hold a weak reference to their broadcaster cannot call it.
  void BroadcastUpdate();

private:
  std::set<Listener *> _listeners;
};

/**
 * The listener for mesh updates. Classes should inherit from this, implement UpdateReceived().
 */
class Listener
{
public:
  /// Initializer.
  Listener(Broadcaster * broadcaster = nullptr) { SetBroadcaster(broadcaster); }

  /// Update the broadcaster. Will replace existing broadcaster.
  void SetBroadcaster(Broadcaster * broadcaster = nullptr) { AddBroadcaster(broadcaster); }

  /// Destructor ensures that we are removed from the broadcaster.
  virtual ~Listener() { RemoveBroadcaster(); }

  /// Override in derived classes.
  virtual void UpdateReceived() { MFEM_ABORT("Not implemented."); }

protected:
  /// Add a new broadcaster.
  void AddBroadcaster(Broadcaster * broadcaster)
  {
    RemoveBroadcaster();

    _broadcaster = broadcaster;
    _broadcaster->AddListener(this);
  }

  /// Remove from broadcasts. Must be called during destruction.
  void RemoveBroadcaster()
  {
    if (_broadcaster)
    {
      _broadcaster->RemoveListener(this);
      _broadcaster = nullptr;
    }
  }

private:
  /// Hold a pointer to our broadcaster. For now, assume that the broadcaster will have
  /// a lifetime at least equal to that of the longest-living listener. We may need to
  /// implement a shared pointer to ensure that the member variable is always valid.
  Broadcaster * _broadcaster{nullptr};
};

}
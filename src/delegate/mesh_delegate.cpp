#include "mesh_delegate.hpp"

namespace hephaestus
{
void
Broadcaster::BroadcastUpdate()
{
  for (const auto & listener : _listeners)
  {
    listener->UpdateReceived();
  }
}
}

#pragma once
#include <set>
#include <memory>

// TODO: - this only allows one type of subscription.

namespace hephaestus
{
/**
 * Classes which inherit from Subscriber can respond to an update from a subscribed class.
 */
class Subscriber
{
public:
  Subscriber() = default;
  ~Subscriber() = default;

  /**
   * Update after receiving some subscribed broadcast.
   */
  virtual void OnDidReceiveNotification() = 0;
};

/**
 * Classes which inherit from Broadcaster can subscribe over classes which inherit from Subscriber.
 * This allows these classes to be notified by an update.
 */
class Broadcaster
{
public:
  Broadcaster() = default;
  ~Broadcaster() { RemoveSubscribers(); }

  /**
   * Notify all subscribers.
   */
  virtual void Broadcast()
  {
    // Call update method for each subscriber.
    for (const auto & subscriber : _subscribers)
    {
      subscriber->OnDidReceiveNotification();
    }
  }

  /**
   * Add a new subscriber.
   */
  void AddSubscriber(std::shared_ptr<Subscriber> subscriber)
  {
    if (_subscribers.count(subscriber) == 0)
    {
      _subscribers.insert(std::move(subscriber));
    }
  }

  /**
   * Remove all subscribers.
   */
  void RemoveSubscribers() { _subscribers.clear(); }

private:
  // TODO: - could this result in strong reference cycles.
  std::set<std::shared_ptr<Subscriber>> _subscribers;
};
}
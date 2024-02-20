#pragma once
#include "timer.hpp"
#include <iostream>

namespace hephaestus
{

class Logger
{

private:
  hephaestus::Timer<std::chrono::seconds> _timer;
  std::string _log_name;
  bool _active = false;

public:
  Logger(const std::string logname, bool active) : _log_name(logname), _active(active){};

  template <typename T>
  friend void operator<<(const Logger & logger, const T stream);
};

template <typename T>
inline void
operator<<(const Logger & logger, const T stream)
{
  if (logger._active)
    std::cout << logger._timer << " : Hephaestus Log " + logger._log_name + " : " << stream << "\n";
}

} // namespace hephaestus

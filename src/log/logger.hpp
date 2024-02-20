#pragma once

#include "timer.hpp"
#include "mfem.hpp"
#include <iostream>

namespace hephaestus
{

class Logger
{

private:
  hephaestus::Timer<std::chrono::seconds> _timer;
  std::vector<hephaestus::Timer<std::chrono::seconds>> _section_timers;
  std::string _log_name;
  std::ostream & _outstream = std::cout;
  bool _active = false;

public:
  Logger(const std::string logname, bool active) : _log_name(logname), _active(active){};

  void LogPush(std::string timername);
  void LogPop();

  template <typename T>
  friend void operator<<(const Logger & logger, const T stream);
};

template <typename T>
inline void
operator<<(const Logger & logger, const T stream)
{
  if (logger._active && mfem::Mpi::WorldRank() == 0)
    logger._outstream << logger._timer << " : Hephaestus Log " + logger._log_name + " : " << stream
                      << "\n";
}

} // namespace hephaestus

// Our selection of logs
extern hephaestus::Logger LogMessage;
extern hephaestus::Logger LogPerformance;
extern hephaestus::Logger LogDebug;

#include "logger.hpp"

namespace hephaestus
{

void
Logger::LogPush(std::string timername)
{
  _section_timers.emplace_back(6, timername);
}

void
Logger::LogPop()
{

  if (_section_timers.size() == 0)
    mfem::mfem_error("LogPop cannot be matched to a previous LogPush");

  std::stringstream stream;
  stream << "Time in " << _section_timers.back().GetName() << " section was "
         << _section_timers.back();
  *this << stream.str();

  _section_timers.pop_back();
}

}

// Instantiating loggers
hephaestus::Logger LogMessage("Message", true);
hephaestus::Logger LogPerformance("Performance", true);
hephaestus::Logger LogDebug("Debug", true);

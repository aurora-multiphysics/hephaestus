#pragma once

#include <ctime>
#include <chrono>

namespace hephaestus
{

class Timer
{
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;

  TimePoint _timerstart;

  Timer() : _timerstart(Clock::now()) {}

  inline double Microsec()
  {
    auto musec = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - _timerstart);
    return musec.count();
  }
};

} // namespace hephaestus
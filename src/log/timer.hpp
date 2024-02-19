#pragma once

#include <ostream>
#include <chrono>

namespace hephaestus
{

class Timer
{
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;

  TimePoint _timerstart;
  int _precision = 6;

  Timer() : _timerstart(Clock::now()) {}
  Timer(int precision) : _timerstart(Clock::now()), _precision(precision) {}

  [[nodiscard]] inline double Seconds() const
  {
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - _timerstart);
    return sec.count();
  }

  friend std::ostream & operator<<(std::ostream & stream, const Timer & timer);
};

inline std::ostream &
operator<<(std::ostream & stream, const Timer & timer)
{
  stream << timer.Seconds() << " s";
  return stream;
}

} // namespace hephaestus
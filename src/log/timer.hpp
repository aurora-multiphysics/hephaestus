#pragma once

#include <chrono>
#include <ostream>
#include <iomanip>

namespace hephaestus
{

template <typename T = std::chrono::microseconds>
class Timer
{
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;

private:
  TimePoint _timerstart;
  int _precision = 6;

public:
  Timer() : _timerstart(Clock::now()) {}
  Timer(int precision) : _timerstart(Clock::now()), _precision(precision) {}

  [[nodiscard]] inline double Microsec() const
  {
    auto usec = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - _timerstart);
    return usec.count();
  }

  friend std::ostream & operator<<(std::ostream & stream, const Timer & timer);
};

inline std::ostream &
operator<<(std::ostream & stream, const Timer<std::chrono::microseconds> & timer)
{
  stream << std::fixed << std::setprecision(timer._precision) << timer.Microsec() << " μs";
  return stream;
}

inline std::ostream &
operator<<(std::ostream & stream, const Timer<std::chrono::milliseconds> & timer)
{
  stream << std::fixed << std::setprecision(timer._precision) << timer.Microsec() / 1e3 << " ms";
  return stream;
}

inline std::ostream &
operator<<(std::ostream & stream, const Timer<std::chrono::seconds> & timer)
{
  stream << std::fixed << std::setprecision(timer._precision) << timer.Microsec() / 1e6 << " s";
  return stream;
}

} // namespace hephaestus
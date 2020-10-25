#include "performance_utils.h"
Timer::Timer() : time_start(Clock::now()) {}
void Timer::start() { time_start = Clock::now(); }
double Timer::elapsed() {
  std::chrono::time_point<Clock> time_end = Clock::now();
  using Ms = std::chrono::nanoseconds;
  return std::chrono::duration_cast<Ms>(time_end - time_start).count();
}

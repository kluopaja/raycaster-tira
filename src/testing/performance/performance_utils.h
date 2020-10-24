#include <chrono>
class Timer {
 public:
   Timer();
   void start();
   double elapsed();
 private:
  using Clock = std::chrono::high_resolution_clock;
  std::chrono::time_point<Clock> time_start;
};

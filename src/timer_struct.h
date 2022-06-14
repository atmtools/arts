#ifndef timer_struct_h
#define timer_struct_h

#include <chrono>

struct Timer {
  bool running{false};
  bool finished{false};
  std::clock_t cputime_start;
  std::chrono::time_point<std::chrono::high_resolution_clock> realtime_start;
  std::clock_t cputime_end;
  std::chrono::time_point<std::chrono::high_resolution_clock> realtime_end;
};

#endif

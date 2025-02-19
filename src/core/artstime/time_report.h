#pragma once

#include <source_location>

#include "artstime.h"

namespace arts {
struct profiler {
  std::string name;
  Time start;

  profiler(std::source_location loc = std::source_location::current())
      : name(loc.function_name()), start(Time{}) {}

  ~profiler();
};

std::string get_report(bool clear=true);
}  // namespace arts

#if ARTS_PROFILING
#define ARTS_TIME_REPORT arts::profiler _arts_prof_var_name_{};
#else
#define ARTS_TIME_REPORT
#endif

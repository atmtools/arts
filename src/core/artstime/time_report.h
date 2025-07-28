#pragma once

#include <source_location>

#include "artstime.h"

namespace arts {
struct profiler {
  std::string name;
  Time start;

  profiler(std::string&& key);
  profiler(std::source_location loc = std::source_location::current());
  
  profiler(const profiler&) = delete;
  profiler(profiler&&) = delete;
  profiler& operator=(const profiler&) = delete;
  profiler& operator=(profiler&&) = delete;

  ~profiler();
};

std::string get_report(Size min_time = 0, bool clear = true);
}  // namespace arts

#if ARTS_PROFILING
#define ARTS_TIME_REPORT arts::profiler _arts_prof_var_name_{};
#define ARTS_NAMED_TIME_REPORT(_name_var_) arts::profiler _arts_named__prof_var_name_{_name_var_};
#else
#define ARTS_TIME_REPORT
#define ARTS_NAMED_TIME_REPORT(_name_var_)
#endif

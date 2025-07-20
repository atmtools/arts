#include "time_report.h"

#include <chrono>
#include <mutex>
#include <unordered_map>

namespace arts {
namespace {
std::string short_name(const std::string& name) {
  return split(split(name, "(")[0], " ").back();
}
}  // namespace

TimeReport profile_report;
std::mutex mprofile_report;

profiler::profiler(std::string&& key) : name(std::move(key)), start(Time{}) {}

profiler::profiler(std::source_location loc)
    : profiler(short_name(loc.function_name())) {}

profiler::~profiler() {
  const Time end{};
  const int core = arts_omp_get_thread_num();

  // Lock might add extra pause inbetween thread calls, but it is safe
  std::scoped_lock lock{mprofile_report};
  profile_report[core][name].emplace_back(start, end);
}

TimeReport get_report(bool clear) {
  std::scoped_lock lock{mprofile_report};

  TimeReport copy = profile_report;

  if (clear)
    for (auto& [key, value] : profile_report) value.clear();
  return copy;
}
}  // namespace arts
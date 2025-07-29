#include "time_report.h"

#include <arts_omp.h>

#include <mutex>
#include <unordered_map>

namespace arts {
namespace {
std::string short_name(const std::string& name) {
  const auto p1 = split(name, "(");

  const auto s1 = std::span{p1}.first(std::max<Size>(p1.size() - 1, 1));
  String s      = join(s1, "(");
  trim(s);

  const auto p2 = split(s, " ");
  const auto s2 = std::span{p2}.last(std::max<Size>(p2.size() - 1, 1));

  return join(s2, " ");
}

TimeReport profile_report;
std::mutex mprofile_report;
}  // namespace

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
#include "time_report.h"

#include <arts_conversions.h>
#include <arts_omp.h>
#include <mystring.h>

#include <map>
#include <mutex>
#include <print>
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

profiler::profiler(std::string&& key)
    : name(std::move(key)), start(std::chrono::system_clock::now()) {}

profiler::profiler(std::source_location loc)
    : profiler(short_name(loc.function_name())) {}

profiler::~profiler() {
  const time_t end{std::chrono::system_clock::now()};
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

void print_report() {
  const auto report = get_report(false);

  // sort
  std::map<std::string, std::vector<StartEnd>> merged;
  for (const auto& [core, core_report] : report) {
    for (const auto& [name, times] : core_report) {
      auto& merged_times = merged[name];
      merged_times.insert(merged_times.end(), times.begin(), times.end());
    }
  }

  std::println("| Method | Average Time | Min Time | Max Time | # Calls |");
  std::println("|---|---|---|---|---|");
  for (const auto& [name, times] : merged) {
    Numeric total = 0.0;
    Numeric max   = std::numeric_limits<Numeric>::lowest();
    Numeric min   = std::numeric_limits<Numeric>::max();
    for (const auto& [start, end] : times) {
      const Numeric t  = std::chrono::duration<Numeric>(end - start).count();
      total           += t;
      max              = std::max(max, t);
      min              = std::min(min, t);
    }
    const Numeric avg = total / static_cast<Numeric>(times.size());

    const auto [ua, a] = Conversion::metric_prefix(avg);
    const auto [um, m] = Conversion::metric_prefix(min);
    const auto [ux, x] = Conversion::metric_prefix(max);

    std::println(
        "| {} | {:.2f} {}s | {:.2f} {}s | {:.2f} {}s | {} |",
        name,
        a,
        ua,
        m,
        um,
        x,
        ux,
        times.size());
  }
}
}  // namespace arts
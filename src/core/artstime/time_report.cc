#include "time_report.h"

#include <sorting.h>
#include <mutex>

namespace {
using Delta = Time::InternalTimeStep;
std::unordered_map<std::string, std::vector<Delta>> profile_report;
std::mutex m;
};  // namespace

namespace arts {
profiler::~profiler() {
  const Time end{};
  const auto diff = end - start;
  std::scoped_lock lock{m};
  profile_report[name].push_back(diff);
}

std::string get_report(bool clear) {
  if (profile_report.empty()) {
#if ARTS_PROFILING
    return "No profiling data available.\n";
#else
    return "Automatic profiling is disabled and no manual profiling performed.\n";
#endif
  }

  std::scoped_lock lock{m};

  std::vector<std::string> vec{};
  std::vector<Delta> total_time{};
  for (const auto& [name, deltas] : profile_report) {
    const Size N = deltas.size();
    if (N == 0) continue;

    const auto sum = std::accumulate(deltas.begin(), deltas.end(), Delta{});
    const auto avg = sum / N;
    const auto max = *stdr::max_element(deltas);
    const auto min = *stdr::min_element(deltas);

    total_time.push_back(sum);
    vec.push_back(std::format("| {0} | {1}μs | {3}μs | {4}μs | {5}μs | {2} |\n",
                              name,
                              std::chrono::duration_cast<std::chrono::microseconds>(avg).count(),
                              N,
                              std::chrono::duration_cast<std::chrono::microseconds>(min).count(),
                              std::chrono::duration_cast<std::chrono::microseconds>(max).count(),
                              std::chrono::duration_cast<std::chrono::microseconds>(sum).count()));
  }

  if (clear) profile_report.clear();

  bubble_sort_by(
      [&](Index i, Index j) { return total_time[i] < total_time[j]; },
      total_time,
      vec);

  return std::format(
      "| Function | Average time | Min time | Max time | Total time | Times called |\n"
      "| -------- | ------------ | -------- | -------- | ---------- | ------------ |\n"
      "{}",
      vec);
}
}  // namespace arts
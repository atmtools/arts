#include "time_report.h"

#include <mutex>

using Delta = Time::InternalTimeStep;
std::unordered_map<std::string, std::vector<Delta>> profile_report;
std::mutex m;

namespace arts {
profiler::~profiler() {
  const Time end{};
  const auto diff = end - start;
  std::scoped_lock lock{m};
  profile_report[name].push_back(diff);
}

std::string get_report(bool clear) {
  if (profile_report.empty()) return "No profiling data available.\n";

  std::scoped_lock lock{m};

  std::string out{
      "| Function-name | Average Time | Times Called | Min Time | Max Time |\n"};
  auto inc = std::back_inserter(out);
  for (const auto& [name, deltas] : profile_report) {
    const Size N = deltas.size();
    if (N == 0) continue;

    const auto sum = std::accumulate(deltas.begin(), deltas.end(), Delta{});
    const auto avg = sum / N;
    const auto max = *stdr::max_element(deltas);
    const auto min = *stdr::min_element(deltas);
    std::format_to(
        inc, "| {0} | {1} | {2} | {3} | {4} |\n", name, avg, N, min, max);
  }

  if (clear) profile_report.clear();

  return out;
}
}  // namespace arts
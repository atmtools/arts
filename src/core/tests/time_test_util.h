#pragma once

#include <algorithm>
#include <chrono>
#include <map>
#include <numeric>
#include <print>
#include <utility>
#include <vector>

std::map<std::string, std::vector<long long>> time_points;

struct test_timer_t {
  using Clock      = std::chrono::high_resolution_clock;
  using time_point = std::chrono::time_point<Clock>;

  std::string name;
  time_point start;
  test_timer_t(std::string n)
      : name(std::move(n)), start(std::chrono::high_resolution_clock::now()) {}

  ~test_timer_t() {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();
    time_points[name].push_back(duration);
  }
};

inline void print_time_points() {
  for (const auto &entry : time_points) {
    const auto &name  = entry.first;
    const auto &times = entry.second;

    auto min = stdr::min_element(times);
    auto max = stdr::max_element(times);
    auto sum = std::accumulate(times.begin(), times.end(), 0LL);

    if (times.size() == 1)
      std::print("{}: {} μs\n", name, *min);
    else
      std::print("{}: min = {} μs, max = {} μs, avg = {} μs, run {} times\n",
                 name,
                 *min,
                 *max,
                 sum / times.size(),
                 times.size());
  }
}

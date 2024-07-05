#pragma once

#include <artstime.h>

#include <chrono>
#include <iomanip>
#include <string_view>

#include "debug.h"

struct Timing {
  std::string_view name;
  Timing(const char* c) : name(c) {}
  TimeStep dt{};
  template <typename Function>
  void operator()(Function&& f) {
    Time start{};
    f();
    Time end{};
    dt = end - start;
  }
};

inline std::ostream& operator<<(std::ostream& os, const Array<Timing>& vt) {
  for (auto& t : vt) {
    if (t.name.contains('\n') or t.name.contains(' ') or t.name.empty())
      throw std::runtime_error(var_string("bad name: \"", t.name, '"'));
    if (t.name not_eq "dummy") {
      os << std::setprecision(15) << t.name << " " << t.dt.count() << '\n';
    }
  }
  return os;
}

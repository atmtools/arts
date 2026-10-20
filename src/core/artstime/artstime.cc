/**
 * @file   artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to time
*/

#include "artstime.h"

#include <debug.h>
#include <enums.h>
#include <fast_float/fast_float.h>

#include <charconv>
#include <chrono>
#include <cstdlib>
#include <format>
#include <istream>
#include <ostream>
#include <print>

Time::Time() : time(std::chrono::system_clock::now()) {}

Time::Time(std::chrono::year_month_day date, InternalTimeStep time_of_day, std::chrono::minutes offset)
    : time(std::chrono::sys_days{date} + time_of_day - offset), utc_offset(offset) {}

std::chrono::year_month_day Time::date() const {
  return std::chrono::year_month_day{std::chrono::floor<std::chrono::days>(time + utc_offset)};
}

std::chrono::hh_mm_ss<Time::InternalTimeStep> Time::time_of_day() const {
  const auto local = time + utc_offset;
  return std::chrono::hh_mm_ss<InternalTimeStep>{local - std::chrono::floor<std::chrono::days>(local)};
}

String Time::str() const {
  String out = std::format("{:%F %T}", time + utc_offset);

  if (utc_offset != std::chrono::minutes{0}) {
    const auto offset  = std::chrono::abs(utc_offset).count();
    out               += std::format("{}{:02}:{:02}", utc_offset.count() < 0 ? '-' : '+', offset / 60, offset % 60);
  }

  return out;
}

TimeStep Time::seconds_into_day() const { return TimeStep{time_of_day().to_duration()}; }

Time::InternalTimeStep Time::EpochTime() const { return time.time_since_epoch(); }

// Operations
Time::InternalTimeStep Time::operator-(const Time& t) const noexcept { return time - t.time; }

// helpers
Numeric Time::Seconds() const { return std::chrono::duration_cast<TimeStep>(time.time_since_epoch()).count(); }

void    Time::Seconds(Numeric x) { time = decltype(time){std::chrono::round<InternalTimeStep>(TimeStep(x))}; }
Numeric Time::PartOfSecond() const { return std::fmod(Seconds(), 1.0); }

// Conversion
Time::operator Numeric() const { return Seconds(); }

Time::Time(const String& t) {
  auto s = std::istringstream(t);
  s >> *this;
}

TimeStep time_stepper_selection(const String& time_step) {
  std::istringstream x(time_step);
  Index              length;
  String             type;
  x >> length >> type;
  tolower(type);

  const TimeStepType t = to<TimeStepType>(type);

  switch (t) {
    case TimeStepType::hour:   return TimeStep(std::chrono::hours(length));
    case TimeStepType::minute: return TimeStep(std::chrono::minutes(length));
    case TimeStepType::second: return TimeStep(std::chrono::seconds(length));
  }
  return {};
}

Time next_even(const Time& t, const TimeStep& dt) {
  auto dt_internal = std::chrono::duration_cast<Time::InternalTimeStep>(dt);

  if (dt_internal.count()) return t + dt_internal - t.EpochTime() % dt_internal;
  return t;
}

ArrayOfIndex time_steps(const ArrayOfTime& times, const TimeStep& DT) {
  Index N = times.size();
  ARTS_USER_ERROR_IF(N < 2, "Can only find time steps for 2-long or longer time grids");

  // algorithm only works with absolute times
  const TimeStep dt = std::chrono::abs(DT);

  ArrayOfIndex time_steps{0};

  if (N > 2) {
    Time tupp = next_even(times[0], dt);
    time_steps.emplace_back(1);

    for (; time_steps.back() < N; time_steps.back()++) {
      if (not(times[time_steps.back()] < tupp)) {
        tupp = next_even(tupp, dt);
        time_steps.emplace_back(time_steps.back());
      }
    }
  } else {
    time_steps.emplace_back(N);
  }
  return time_steps;
}

std::ostream& operator<<(std::ostream& os, const Time& t) { return os << t.str(); }

namespace {
template <typename T> T parse_number(std::string_view s, std::string_view what, std::string_view text) {
  T x{};

  const auto [ptr, ec] = [&]() {
    if constexpr (std::floating_point<T>)
      return fast_float::from_chars(s.data(), s.data() + s.size(), x);
    else
      return std::from_chars(s.data(), s.data() + s.size(), x);
  }();

  ARTS_USER_ERROR_IF(ec != std::errc{} or ptr != s.data() + s.size() or s.empty(),
                     "Cannot read the {} from \"{}\" in the time \"{}\"",
                     what,
                     s,
                     text)

  return x;
}
}  // namespace

std::istream& operator>>(std::istream& is, Time& t) {
  using namespace std::chrono;

  String ymd_str, hms_str;
  is >> ymd_str >> hms_str;

  const String text = ymd_str + ' ' + hms_str;

  const auto YMD = split(ymd_str, "-");
  ARTS_USER_ERROR_IF(YMD.size() != 3,
                     R"(The time "{}" is not "YYYY-MM-DD hh:mm:ss", optionally with fractional seconds and an offset
from UTC, e.g., "2000-01-01 12:00:00.5+09:00".  Times without an offset are UTC.)",
                     text)

  // The offset from UTC, "Z", "+hh:mm" or "-hh:mm", at the end of the time of day
  std::string_view hms_view{hms_str};
  minutes          offset{0};
  if (hms_view.ends_with('Z')) {
    hms_view.remove_suffix(1);
  } else if (const auto i = hms_view.find_first_of("+-"); i != std::string_view::npos) {
    const std::string_view off = hms_view.substr(i + 1);
    const auto             c   = off.find(':');
    ARTS_USER_ERROR_IF(
        c == std::string_view::npos, "The offset from UTC must be \"+hh:mm\" or \"-hh:mm\" in \"{}\"", text)

    offset = hours{parse_number<int>(off.substr(0, c), "offset hours", text)} +
             minutes{parse_number<int>(off.substr(c + 1), "offset minutes", text)};
    if (hms_view[i] == '-') offset = -offset;
    hms_view = hms_view.substr(0, i);
  }

  const auto HMS = split(String{hms_view}, ":");
  ARTS_USER_ERROR_IF(HMS.size() != 3,
                     R"(The time "{}" is not "YYYY-MM-DD hh:mm:ss", optionally with fractional seconds and an offset
from UTC, e.g., "2000-01-01 12:00:00.5+09:00".  Times without an offset are UTC.)",
                     text)

  const year_month_day date{year{parse_number<int>(YMD[0], "year", text)},
                            month{parse_number<unsigned>(YMD[1], "month", text)},
                            day{parse_number<unsigned>(YMD[2], "day", text)}};
  const int            hour   = parse_number<int>(HMS[0], "hour", text);
  const int            minute = parse_number<int>(HMS[1], "minute", text);
  const Numeric        second = parse_number<Numeric>(HMS[2], "second", text);

  ARTS_USER_ERROR_IF(not date.ok(), "The date of the time \"{}\" does not exist", text)
  ARTS_USER_ERROR_IF(hour < 0 or hour > 23 or minute < 0 or minute > 59 or second < 0 or second >= 60,
                     "The time of day of the time \"{}\" is not within 00:00:00 to 23:59:59",
                     text)
  ARTS_USER_ERROR_IF(abs(offset) >= days{1}, "The offset from UTC of the time \"{}\" is not less than a day", text)

  t = Time{date, hours{hour} + minutes{minute} + round<Time::InternalTimeStep>(TimeStep{second}), offset};

  return is;
}

Time mean_time(const ArrayOfTime& ts, Index s, Index E) {
  Index e = 0;
  if (e == -1) e = ts.size();
  ARTS_USER_ERROR_IF(
      e < 0 or e > static_cast<Index>(ts.size()), "Bad last index, valid options are [-1, ts.size()], got: {}", E);

  ARTS_USER_ERROR_IF(
      s < 0 or s > static_cast<Index>(ts.size()), "Bad first index, valid options are [0, ts.size()], got: {}", s);

  Time::InternalTimeStep dt(0);
  for (Index i = s + 1; i < e; i++) dt += (ts[i] - ts[s]) / (e - s);
  return ts[s] + dt;
}

TimeStep median(ArrayOfTimeStep dt) {
  const auto n = dt.size();
  if (n) {
    std::sort(dt.begin(), dt.end());
    if (n % 2) return dt[n / 2];
    return (dt[(n - 1) / 2] + dt[n / 2]) / 2;
  }
  return TimeStep(0);
}

TimeStep mean(const ArrayOfTimeStep& dt) {
  TimeStep   t(0);
  const auto n = dt.size();
  for (std::size_t i = 0; i < n; i++) { t += dt[i] / n; }
  return t;
}

std::ostream& operator<<(std::ostream& os, const TimeStep& dt) { return os << dt.count() << " seconds"; }

std::ostream& operator<<(std::ostream& os, const ArrayOfTime& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfTime& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfTimeStep& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

DebugTime::~DebugTime() {
#pragma omp critical
  std::print(stderr, "{}: {}\n", msg, Time{} - start);
}

void xml_io_stream<Time>::write(std::ostream& os, const Time& x, bofstream*, std::string_view name) {
  XMLTag tag(type_name, "name", name, "version", x.Version());
  tag.write_to_stream(os);
  std::println(os, "{:IO}", x);
  tag.write_to_end_stream(os);
}

void xml_io_stream<Time>::read(std::istream& is, Time& x, bifstream*) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  is >> x;

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
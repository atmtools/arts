/**
 * @file   artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to time
*/

#include "artstime.h"

#include <debug.h>
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
  const auto offset = std::chrono::abs(utc_offset).count();
  return std::format(
      "{:%F %T}{}{:02}:{:02}", time + utc_offset, utc_offset.count() < 0 ? '-' : '+', offset / 60, offset % 60);
}

Time::InternalTimeStep Time::EpochTime() const { return time.time_since_epoch(); }

// Operations
Time::InternalTimeStep Time::operator-(const Time& t) const noexcept { return time - t.time; }

Numeric Time::Seconds() const { return std::chrono::duration_cast<TimeStep>(time.time_since_epoch()).count(); }

void Time::Seconds(Numeric x) { time = decltype(time){std::chrono::round<InternalTimeStep>(TimeStep(x))}; }

Time::Time(const String& t) {
  auto s = std::istringstream(t);
  s >> *this;
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

  if (ec != std::errc{} or ptr != s.data() + s.size() or s.empty())
    arts::user_error("Cannot read the {} from \"{}\" in the time \"{}\"", what, s, text);

  return x;
}
}  // namespace

std::istream& operator>>(std::istream& is, Time& t) {
  using namespace std::chrono;

  // The date and the time of day are separated by a space or, as in RFC 3339, by "T"
  String ymd_str, hms_str;
  char   separator = ' ';
  is >> ymd_str;
  if (const auto i = ymd_str.find_first_of("Tt"); i != String::npos) {
    separator = ymd_str[i];
    hms_str   = ymd_str.substr(i + 1);
    ymd_str.resize(i);
  } else {
    is >> hms_str;
  }

  const String text = ymd_str + separator + hms_str;

  // The offset from UTC, "Z", "+hh:mm" or "-hh:mm", at the end of the time of day
  std::string_view hms_view{hms_str};
  minutes          offset{0};
  if (hms_view.ends_with('Z') or hms_view.ends_with('z')) {
    hms_view.remove_suffix(1);
  } else if (const auto i = hms_view.find_first_of("+-"); i != std::string_view::npos) {
    const std::string_view off = hms_view.substr(i + 1);
    const auto             c   = off.find(':');
    if (c == std::string_view::npos)
      arts::user_error("The offset from UTC must be \"+hh:mm\" or \"-hh:mm\" in \"{}\"", text);

    offset = hours{parse_number<int>(off.substr(0, c), "offset hours", text)} +
             minutes{parse_number<int>(off.substr(c + 1), "offset minutes", text)};
    if (hms_view[i] == '-') offset = -offset;
    hms_view = hms_view.substr(0, i);
  }

  const auto YMD = split(ymd_str, "-");
  const auto HMS = split(String{hms_view}, ":");
  if (YMD.size() != 3 or HMS.size() != 3)
    arts::user_error(R"(The time "{}" is not "YYYY-MM-DD hh:mm:ss" or "YYYY-MM-DDThh:mm:ss", optionally with
fractional seconds and an offset from UTC, "Z", "+hh:mm" or "-hh:mm", e.g., "2000-01-01T12:00:00.5+09:00".
Times without an offset are UTC.)",
                     text);

  const year_month_day date{year{parse_number<int>(YMD[0], "year", text)},
                            month{parse_number<unsigned>(YMD[1], "month", text)},
                            day{parse_number<unsigned>(YMD[2], "day", text)}};
  const int            hour   = parse_number<int>(HMS[0], "hour", text);
  const int            minute = parse_number<int>(HMS[1], "minute", text);
  const Numeric        second = parse_number<Numeric>(HMS[2], "second", text);

  if (not date.ok()) arts::user_error("The date of the time \"{}\" does not exist", text);
  if (hour < 0 or hour > 23 or minute < 0 or minute > 59 or second < 0 or second >= 60)
    arts::user_error("The time of day of the time \"{}\" is not within 00:00:00 to 23:59:59", text);
  if (abs(offset) >= days{1}) arts::user_error("The offset from UTC of the time \"{}\" is not less than a day", text);

  t = Time{date, hours{hour} + minutes{minute} + round<Time::InternalTimeStep>(TimeStep{second}), offset};

  return is;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfTime& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfTime& a) {
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
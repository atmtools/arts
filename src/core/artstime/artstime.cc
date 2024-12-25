/**
 * @file   artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to time
*/

#include "artstime.h"

#include <enumsTimeStepType.h>

#include <charconv>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <print>

#include "debug.h"

Time::Time() : time(std::chrono::system_clock::now()) {}

Time::Time(std::time_t t) : time(std::chrono::system_clock::from_time_t(t)) {}

Time::Time(std::tm t) : Time(std::mktime(&t)) {}

// Conversions
std::time_t Time::toTimeT() const {
  return std::chrono::system_clock::to_time_t(time);
}

std::tm Time::toStruct() const {
  std::time_t x = toTimeT();
  std::tm y;
#ifdef _MSC_VER
  errno_t err = localtime_s(&y, &x);
  ARTS_USER_ERROR_IF(err != 0, "Cannot construct time struct")
#else
  tm* z = localtime_r(&x, &y);
  ARTS_USER_ERROR_IF(not z, "Cannot construct time struct")
#endif
  return y;
}

std::tm Time::toGMTStruct() const {
  std::time_t x = toTimeT();
  std::tm y;
#ifdef _MSC_VER
  errno_t err = gmtime_s(&y, &x);
  ARTS_USER_ERROR_IF(err != 0, "Cannot construct time struct")
#else
  tm* z = gmtime_r(&x, &y);
  ARTS_USER_ERROR_IF(not z, "Cannot construct time struct")
#endif
  return y;
}

TimeStep Time::seconds_into_day() const {
  std::tm x = toStruct();
  return TimeStep(x.tm_hour * 3600 + x.tm_min * 60 + x.tm_sec + PartOfSecond());
}

Time::InternalTimeStep Time::EpochTime() const {
  return time.time_since_epoch();
}

// Operations
Time::InternalTimeStep Time::operator-(const Time& t) const noexcept {
  return time - t.time;
}

bool Time::operator<(const Time& t) const noexcept { return time < t.time; }

bool Time::operator==(const Time& t) const noexcept { return time == t.time; }

bool Time::operator!=(const Time& t) const noexcept {
  return not this->operator==(t);
}

bool Time::operator<=(const Time& t) const noexcept {
  return this->operator<(t) or this->operator==(t);
}

bool Time::operator>(const Time& t) const noexcept {
  return not this->operator<=(t);
}

bool Time::operator>=(const Time& t) const noexcept {
  return this->operator>(t) or this->operator==(t);
}

// helpers
Numeric Time::Seconds() const {
  return std::chrono::duration_cast<TimeStep>(time.time_since_epoch()).count();
}

void Time::Seconds(Numeric x) { operator+=(TimeStep(x - Seconds())); }
Numeric Time::PartOfSecond() const { return std::fmod(Seconds(), 1.0); }

// Conversion
Time::operator Numeric() const { return Seconds(); }

Time::Time(const String& t) {
  auto s = std::istringstream(t);
  s >> *this;
}

TimeStep time_stepper_selection(const String& time_step) {
  std::istringstream x(time_step);
  Index length;
  String type;
  x >> length >> type;
  tolower(type);

  const TimeStepType t = to<TimeStepType>(type);

  switch (t) {
    case TimeStepType::hour:
      return TimeStep(std::chrono::hours(length));
    case TimeStepType::minute:
      return TimeStep(std::chrono::minutes(length));
    case TimeStepType::second:
      return TimeStep(std::chrono::seconds(length));
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
  ARTS_USER_ERROR_IF(
      N < 2, "Can only find time steps for 2-long or longer time grids");

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

std::ostream& operator<<(std::ostream& os, const Time& t) {
  // FIXME: C++20 has much better calendar handling (BUT NOT YET...)
  std::tm x = t.toStruct();

  // Deal with seconds
  std::array<char, 2 + 1 + 9 + 100> sec;
  Numeric seconds = Numeric(x.tm_sec) + t.PartOfSecond();
  snprintf(sec.data(), sec.size(), "%.9lf", seconds);

  // Print based on std::tm specs
  return os << 1900 + x.tm_year << '-' << std::setfill('0') << std::setw(2)
            << 1 + x.tm_mon << '-' << std::setfill('0') << std::setw(2)
            << x.tm_mday << ' ' << std::setfill('0') << std::setw(2)
            << x.tm_hour << ':' << std::setfill('0') << std::setw(2) << x.tm_min
            << ':' << std::setfill('0') << std::setw(12) << sec.data()
            << std::setfill(' ') << std::setw(0);
}

std::istream& operator>>(std::istream& is, Time& t) {
  String ymd, hms;
  is >> ymd >> hms;

  const auto YMD = split(ymd, "-");
  const auto HMS = split(hms, ":");

  ARTS_USER_ERROR_IF(
      YMD.size() not_eq 3 or HMS.size() not_eq 3,
      "Time stream must look like \"year-month-day hour:min:seconds\"\n"
      "\"year-month-day\"   looks like: \"{}\"\n"
      "\"hour:min:seconds\" looks like: \"{}\"",
      ymd,
      hms);

  // FIXME: C++20 has much better calendar (BUT NOT YET...)
  int year{}, month{}, day{};
  auto res_year =
      std::from_chars(YMD[0].c_str(), YMD[0].c_str() + YMD[0].size(), year);
  auto res_mon =
      std::from_chars(YMD[1].c_str(), YMD[1].c_str() + YMD[1].size(), month);
  auto res_day =
      std::from_chars(YMD[2].c_str(), YMD[2].c_str() + YMD[2].size(), day);

  ARTS_USER_ERROR_IF(std::make_error_code(res_year.ec) or
                         std::make_error_code(res_mon.ec) or
                         std::make_error_code(res_day.ec),
                     "Cannot understand time point for year-month-day: \"{}\"",
                     ymd)
  ARTS_USER_ERROR_IF(year < 1900,
                     "We cannot yet support times before the year 1900")

  int hour{}, minute{};
  Numeric sec{};
  auto res_hour =
      std::from_chars(HMS[0].c_str(), HMS[0].c_str() + HMS[0].size(), hour);
  auto res_min =
      std::from_chars(HMS[1].c_str(), HMS[1].c_str() + HMS[1].size(), minute);
  auto res_sec = fast_float::from_chars(
      HMS[2].c_str(), HMS[2].c_str() + HMS[2].size(), sec);

  ARTS_USER_ERROR_IF(
      std::make_error_code(res_hour.ec) or std::make_error_code(res_min.ec) or
          std::make_error_code(res_sec.ec),
      "Cannot understand time point for hour:minute:second in: \"{}\"",
      hms)

  std::tm tm_struct{};
  tm_struct.tm_year  = year - 1900;
  tm_struct.tm_mon   = month - 1;
  tm_struct.tm_mday  = day;
  tm_struct.tm_hour  = hour;
  tm_struct.tm_hour  = hour;
  tm_struct.tm_min   = minute;
  tm_struct.tm_isdst = -1;
  t                  = Time{tm_struct} + TimeStep{sec};

  return is;
}

Time mean_time(const ArrayOfTime& ts, Index s, Index E) {
  Index e = 0;
  if (e == -1) e = ts.size();
  ARTS_USER_ERROR_IF(
      e < 0 or e > static_cast<Index>(ts.size()),
      "Bad last index, valid options are [-1, ts.size()], got: {}",
      E);

  ARTS_USER_ERROR_IF(
      s < 0 or s > static_cast<Index>(ts.size()),
      "Bad first index, valid options are [0, ts.size()], got: {}",
      s);

  Time::InternalTimeStep dt(0);
  for (Index i = s + 1; i < e; i++) dt += (ts[i] - ts[s]) / (e - s);
  return ts[s] + dt;
}

Vector time_vector(const ArrayOfTime& times) {
  Vector t(times.size());
  for (Size i = 0; i < times.size(); i++) t[i] = Numeric(times[i]);
  return t;
}

ArrayOfTime time_vector(const Vector& times) {
  ArrayOfTime t(times.size());
  for (Size i = 0; i < times.size(); i++) t[i].Seconds(times[i]);
  return t;
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
  TimeStep t(0);
  const auto n = dt.size();
  for (std::size_t i = 0; i < n; i++) {
    t += dt[i] / n;
  }
  return t;
}

std::ostream& operator<<(std::ostream& os, const TimeStep& dt) {
  return os << dt.count() << " seconds";
}

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
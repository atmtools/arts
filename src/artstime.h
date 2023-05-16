/**
 * @file   artstime.h
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to time in ARTS
*/

#ifndef ARTSTIME_H
#define ARTSTIME_H

#include <chrono>
#include <cmath>
#include <string_view>

#include "array.h"
#include "debug.h"
#include "matpack_data.h"
#include "mystring.h"

/** A duration of time, 1 full tick should be 1 second */
using TimeStep = std::chrono::duration<Numeric>;

/** Class to handle time in ARTS */
struct Time {
  std::chrono::system_clock::time_point time;
  using InternalTimeStep = decltype(time)::duration;

  // Version will be updated when C++20 datetime is available... Version 1 will still assume local time at that time
  [[nodiscard]] Index Version() const noexcept { return 1; }

  // Construction
  Time() : time(std::chrono::system_clock::now()) {}
  explicit Time(std::time_t t)
      : time(std::chrono::system_clock::from_time_t(t)) {}
  explicit Time(std::tm t) : Time(std::mktime(&t)) {}
  explicit Time(const String& t);

  // Conversions
  [[nodiscard]] std::time_t toTimeT() const {
    return std::chrono::system_clock::to_time_t(time);
  }
  [[nodiscard]] std::tm toStruct() const {
    std::time_t x = toTimeT();
    std::tm y;
    tm* z = localtime_r(&x, &y);
    ARTS_USER_ERROR_IF(not z, "Cannot construct time struct")
    return y;
  }
  [[nodiscard]] std::tm toGMTStruct() const {
    std::time_t x = toTimeT();
    std::tm y;
    tm* z = gmtime_r(&x, &y);
    ARTS_USER_ERROR_IF(not z, "Cannot construct time struct")
    return y;
  }
  [[nodiscard]] TimeStep seconds_into_day() const {
    std::tm x = toStruct();
    return TimeStep(x.tm_hour * 3600 + x.tm_min * 60 + x.tm_sec +
                    PartOfSecond());
  }
  [[nodiscard]] InternalTimeStep EpochTime() const {
    return time.time_since_epoch();
  }

  // Operations
  InternalTimeStep operator-(const Time& t) const noexcept {
    return time - t.time;
  }
  bool operator<(const Time& t) const noexcept { return time < t.time; }
  bool operator==(const Time& t) const noexcept { return time == t.time; }
  bool operator!=(const Time& t) const noexcept {
    return not this->operator==(t);
  }
  bool operator<=(const Time& t) const noexcept {
    return this->operator<(t) or this->operator==(t);
  }
  bool operator>(const Time& t) const noexcept {
    return not this->operator<=(t);
  }
  bool operator>=(const Time& t) const noexcept {
    return this->operator>(t) or this->operator==(t);
  }
  template <typename T, typename R>
  Time& operator+=(const std::chrono::duration<T, R>& dt) {
    time += std::chrono::duration_cast<InternalTimeStep>(dt);
    return *this;
  }
  template <typename T, typename R>
  Time& operator-=(const std::chrono::duration<T, R>& dt) {
    time -= std::chrono::duration_cast<InternalTimeStep>(dt);
    return *this;
  }
  template <typename T, typename R>
  Time operator+(const std::chrono::duration<T, R>& dt) const {
    return (Time(*this) += dt);
  }
  template <typename T, typename R>
  Time operator-(const std::chrono::duration<T, R>& dt) const {
    return (Time(*this) -= dt);
  }

  // helpers
  [[nodiscard]] Numeric Seconds() const {
    return std::chrono::duration_cast<TimeStep>(time.time_since_epoch())
        .count();
  }
  void Seconds(Numeric x) { operator+=(TimeStep(x - Seconds())); }
  [[nodiscard]] Numeric PartOfSecond() const {
    return std::fmod(Seconds(), 1.0);
  }

  // Conversion
  explicit operator Numeric() const { return Seconds(); }

  friend std::ostream& operator<<(std::ostream& os, const Time& t);

  friend std::istream& operator>>(std::istream& is, Time& t);
};  // Time

/** List of times */
using ArrayOfTime = Array<Time>;

/** List of times */
using ArrayOfArrayOfTime = Array<ArrayOfTime>;

/** List of time steps */
using ArrayOfTimeStep = Array<TimeStep>;

/** Debug output for duration */
std::ostream& operator<<(std::ostream& os, const TimeStep& dt);

/** Returns a time step from valid string
 * 
 * The string should look like: "X type".
 * where X is parsed as a double and where
 * valid string type(s) are:
 * "hour", "hours", "h": Returns a X hour time step
 * "minute", "minutes", "min": Returns a X minute time step
 * "second", "seconds", "s": Returns a X seconds time step
 * 
 * @param[in] time_step A time step of format "X type"
 * @return The time step as TimeStep
 */
TimeStep time_stepper_selection(const String& time_step);

/** Returns the next time after t with an even time-step
 * 
 * @param[in] t A time
 * @param[in] dt A duration of time
 * @return Next even time, e.g., 14:14:00 with dt as 10 minutes gives 14:20:00
 */
Time next_even(const Time& t, const TimeStep& dt);

/** Finds the index matching demands in a list of times
 * 
 * The first index is 0 and the second index is the start of
 * the first even period of the given stepsize
 * 
 * The last index is times.nelem().  If output has 1 element,
 * no range was found matching the criteria.
 * 
 * @param[in] times Times sorted in ascending order
 * @param[in] dt A duration of time
 * @return Starting index of the time-series
 */
ArrayOfIndex time_steps(const ArrayOfTime& times, const TimeStep& dt);

/** Computes the average time in a list
 * 
 * @param[in] ts A list of time
 * @param[in] s A starting index; valid range [0, ts.nelem())
 * @param[in] e The end+1 index; valid range [-1, ts.nelem()]; -1 is treated as ts.nelem()
 */
Time mean_time(const ArrayOfTime& ts, Index s = 0, Index e = -1);

/** Converts from each Time to seconds and returns as Vector
 * 
 * @param[in] times Times
 * @return Vector of Time->Seconds() calls
 */
Vector time_vector(const ArrayOfTime& times);

/** Converts from each Vector entry by seconds and returns as ArrayOfTime
 * 
 * @param[in] times Times
 * @return ArrayOfTime from seconds
 */
ArrayOfTime time_vector(const Vector& times);

/** Returns the median time step
 * 
 * Takes vector by copy to sort
 * 
 * @param[in] times Times
 * @return median time step
 */
TimeStep median(ArrayOfTimeStep);

/** Returns the mean time step
 * 
 * @param[in] times Times
 * @return mean time step
 */
TimeStep mean(const ArrayOfTimeStep&);

//! Used to debug execution time, prints msg+time on destruction to std::cerr
struct DebugTime {
  Time start{};
  std::string_view msg;
  DebugTime(const std::string_view s="Time") : msg(s) {}
  ~DebugTime() {
    #pragma omp critical
    std::cerr << msg << ':' << ' ' << Time{} - start << '\n';
  }
};

#endif  // ARTSTIME_H

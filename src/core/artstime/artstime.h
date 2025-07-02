#pragma once

#include <matpack.h>
#include <mystring.h>
#include <xml.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <string_view>

/** A duration of time, 1 full tick should be 1 second */
using TimeStep = std::chrono::duration<Numeric>;

/** Class to handle time in ARTS */
struct Time {
  std::chrono::time_point<std::chrono::system_clock> time;
  using InternalTimeStep = decltype(time)::duration;

  // Version will be updated when C++20 datetime is available... Version 1 will still assume local time at that time
  [[nodiscard]] static constexpr Index Version() noexcept { return 1; }

  // Construction
  Time();
  explicit Time(std::time_t t);
  explicit Time(std::tm t);
  explicit Time(const String& t);

  // Conversions
  [[nodiscard]] std::time_t toTimeT() const;
  [[nodiscard]] std::tm toStruct() const;
  [[nodiscard]] std::tm toGMTStruct() const;

  [[nodiscard]] TimeStep seconds_into_day() const;

  [[nodiscard]] InternalTimeStep EpochTime() const;

  // Operations
  InternalTimeStep operator-(const Time& t) const noexcept;
  bool operator<(const Time& t) const noexcept;
  bool operator==(const Time& t) const noexcept;
  bool operator!=(const Time& t) const noexcept;
  bool operator<=(const Time& t) const noexcept;
  bool operator>(const Time& t) const noexcept;
  bool operator>=(const Time& t) const noexcept;
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
  [[nodiscard]] Numeric Seconds() const;
  void Seconds(Numeric x);
  [[nodiscard]] Numeric PartOfSecond() const;

  // Conversion
  explicit operator Numeric() const;

  friend std::ostream& operator<<(std::ostream& os, const Time& t);

  friend std::istream& operator>>(std::istream& is, Time& t);
};  // Time

/** Debug output for duration */
std::ostream& operator<<(std::ostream& os, const TimeStep& dt);

/** List of times */
using ArrayOfTime = Array<Time>;

/** List of times */
using ArrayOfArrayOfTime = Array<ArrayOfTime>;

/** List of time steps */
using ArrayOfTimeStep = Array<TimeStep>;

std::ostream& operator<<(std::ostream& os, const ArrayOfTime& a);

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfTime& a);

std::ostream& operator<<(std::ostream& os, const ArrayOfTimeStep& a);

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
 * The last index is times.size().  If output has 1 element,
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
 * @param[in] s A starting index; valid range [0, ts.size())
 * @param[in] e The end+1 index; valid range [-1, ts.size()]; -1 is treated as ts.size()
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
  DebugTime(const std::string_view s) : msg(s) {}
  DebugTime(const char* s = "Time") : msg(s) {}
  DebugTime(std::string&&) =
      delete;  // Class keeps string-view, cannot move from a string
  ~DebugTime();
};

template <>
struct std::formatter<Time> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Time& v, FmtContext& ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(ctx.out(), "{}{}{}", quote, v.time, quote);
  }
};

template <>
struct xml_io_stream<Time> {
  static constexpr std::string_view type_name = "Time"sv;

  static void write(std::ostream& os,
                    const Time& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Time& x,
                   bifstream* pbifs = nullptr);
};

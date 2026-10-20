#pragma once

#include <mystring.h>
#include <xml.h>

#include <chrono>
#include <string_view>

/** A duration of time, 1 full tick should be 1 second */
using TimeStep = std::chrono::duration<Numeric>;

/** A point in time, and the timezone that its calendar date and time of day is in
 *
 * The timezone is a fixed offset from UTC, by default UTC itself.  Times are
 * read and printed as "YYYY-MM-DD hh:mm:ss", optionally with fractional
 * seconds and the offset, e.g., "2000-01-01 12:00:00.5+09:00".  Times without
 * an offset are UTC.  Times compare as points in time, whatever their timezones.
 */
struct Time {
  std::chrono::time_point<std::chrono::system_clock> time;
  using InternalTimeStep = decltype(time)::duration;

  //! The timezone as an offset from UTC
  std::chrono::minutes utc_offset{0};

  [[nodiscard]] static constexpr Index Version() noexcept { return 2; }

  // Construction
  Time();
  explicit Time(const String& t);

  //! The time at a calendar date and time of day in the timezone utc_offset
  Time(std::chrono::year_month_day date, InternalTimeStep time_of_day, std::chrono::minutes utc_offset = {});

  //! The calendar date in the timezone of the time
  [[nodiscard]] std::chrono::year_month_day date() const;

  //! The time of day in the timezone of the time
  [[nodiscard]] std::chrono::hh_mm_ss<InternalTimeStep> time_of_day() const;

  //! The date and time in the timezone of the time, with its offset unless UTC, e.g., "2000-01-01 12:00:00.000000+09:00"
  [[nodiscard]] String str() const;

  [[nodiscard]] TimeStep seconds_into_day() const;

  [[nodiscard]] InternalTimeStep EpochTime() const;

  // Operations
  InternalTimeStep operator-(const Time& t) const noexcept;

  [[nodiscard]] bool operator==(const Time& t) const noexcept { return time == t.time; }
  [[nodiscard]] auto operator<=>(const Time& t) const noexcept { return time <=> t.time; }

  template <typename T, typename R> Time& operator+=(const std::chrono::duration<T, R>& dt) {
    time += std::chrono::duration_cast<InternalTimeStep>(dt);
    return *this;
  }
  template <typename T, typename R> Time& operator-=(const std::chrono::duration<T, R>& dt) {
    time -= std::chrono::duration_cast<InternalTimeStep>(dt);
    return *this;
  }
  template <typename T, typename R> Time operator+(const std::chrono::duration<T, R>& dt) const {
    return (Time(*this) += dt);
  }
  template <typename T, typename R> Time operator-(const std::chrono::duration<T, R>& dt) const {
    return (Time(*this) -= dt);
  }

  // helpers
  [[nodiscard]] Numeric Seconds() const;
  void                  Seconds(Numeric x);
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
  Time             start{};
  std::string_view msg;
  DebugTime(const std::string_view s) : msg(s) {}
  DebugTime(const char* s = "Time") : msg(s) {}
  DebugTime(std::string&&) = delete;  // Class keeps string-view, cannot move from a string
  ~DebugTime();
};

template <> struct std::formatter<Time> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const Time& v, FmtContext& ctx) const {
    const std::string_view quote = tags.quote();
    return tags.format(ctx, quote, v.str(), quote);
  }
};

template <> struct xml_io_stream<Time> {
  static constexpr std::string_view type_name = "Time"sv;

  static void write(std::ostream& os, const Time& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, Time& x, bifstream* pbifs = nullptr);
};

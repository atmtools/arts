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
 * read as "YYYY-MM-DD hh:mm:ss" or, as in RFC 3339, "YYYY-MM-DDThh:mm:ss",
 * optionally with fractional seconds and the offset, "Z", "+hh:mm" or "-hh:mm",
 * e.g., "2000-01-01T12:00:00.5+09:00", where times without an offset are UTC.
 * Times compare as points in time, whatever their timezones.
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

  //! The date and time in the timezone of the time, with its offset, e.g., "2000-01-01 12:00:00.000000+09:00"
  [[nodiscard]] String str() const;

  [[nodiscard]] InternalTimeStep EpochTime() const;

  // Operations
  InternalTimeStep operator-(const Time& t) const noexcept;

  [[nodiscard]] bool operator==(const Time& t) const noexcept { return time == t.time; }
  [[nodiscard]] auto operator<=>(const Time& t) const noexcept { return time <=> t.time; }

  //! The seconds since the start of Unix time
  [[nodiscard]] Numeric Seconds() const;
  void                  Seconds(Numeric x);

  friend std::ostream& operator<<(std::ostream& os, const Time& t);

  friend std::istream& operator>>(std::istream& is, Time& t);
};  // Time

/** List of times */
using ArrayOfTime = Array<Time>;

/** List of times */
using ArrayOfArrayOfTime = Array<ArrayOfTime>;

std::ostream& operator<<(std::ostream& os, const ArrayOfTime& a);

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfTime& a);

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

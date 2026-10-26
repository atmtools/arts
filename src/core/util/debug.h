/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#include <format_tags.h>

#include <cstdint>
#include <format>
#include <source_location>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <vector>

struct src_location {
  // Snapshot at compile time. Keeping the implementation pointer from
  // std::source_location can expose mis-coalesced local records on GCC/Darwin
  // when a translation unit also emits weak template-static data.
  consteval src_location(std::source_location loc = std::source_location::current())
      : file_name_(loc.file_name()), function_name_(loc.function_name()), line_(loc.line()), column_(loc.column()) {}

  std::string get() const;
  std::string getfunc() const;

 private:
  const char*         file_name_;
  const char*         function_name_;
  std::uint_least32_t line_;
  std::uint_least32_t column_;
};

namespace arts {
std::runtime_error catch_errors(std::logic_error& e, const std::string_view context);
std::runtime_error catch_errors(std::exception& e, const std::string_view context);

//! Whether the format string refers to every one of the n arguments
consteval bool uses_all_arguments(const std::string_view fmt, const std::size_t n) {
  std::vector<bool> used(n, false);
  std::size_t       next = 0;

  const auto field = [&](const std::string_view id) {
    std::size_t i = 0;
    if (id.empty()) {
      i = next++;
    } else {
      for (const char c : id) i = 10 * i + static_cast<std::size_t>(c - '0');
    }
    if (i < n) used[i] = true;
  };

  for (std::size_t i = 0; i < fmt.size(); ++i) {
    if (fmt[i] == '}') {
      ++i;  // "}}", a single '}' is rejected by std::format_string
    } else if (fmt[i] == '{') {
      if (fmt[i + 1] == '{') {
        ++i;
        continue;
      }

      std::size_t j = i + 1;
      while (fmt[j] >= '0' and fmt[j] <= '9') ++j;
      field(fmt.substr(i + 1, j - i - 1));

      // The format specification may hold replacement fields for width and precision
      for (; fmt[j] != '}'; ++j) {
        if (fmt[j] == '{') {
          const std::size_t k = fmt.find('}', j);
          field(fmt.substr(j + 1, k - j - 1));
          j = k;
        }
      }
      i = j;
    }
  }

  for (const bool u : used) {
    if (not u) return false;
  }
  return true;
}

//! Not constexpr: calling it from a consteval function makes the call fail to compile
void the_format_string_does_not_use_every_argument();

//! A format string that remembers where it was written, and that uses all of its arguments
template <typename... Args> struct located_format {
  std::format_string<Args...> fmt;
  src_location                location;

  template <typename S> consteval located_format(const S& s, const src_location loc = std::source_location::current())
      : fmt(s), location(loc) {
    if (not uses_all_arguments(fmt.get(), sizeof...(Args))) the_format_string_does_not_use_every_argument();
  }
};

//! Whether user errors report the C++ source location they were thrown from
bool show_cpp_error_location();

//! Sets whether user errors report the C++ source location they were thrown from
void show_cpp_error_location(bool show);

//! Throws a user error with the message
[[noreturn]] void user_error(const std::string_view msg,
                             const src_location&    location = std::source_location::current());

//! Throws a user error with the formatted message
template <std::formattable<char>... Args>
[[noreturn]] void user_error(const located_format<std::type_identity_t<Args>...> fmt, const Args&... args)
    requires(sizeof...(Args) > 0) {
  user_error(std::vformat(fmt.fmt.get(), std::make_format_args(args...)), fmt.location);
}
}  // namespace arts

#define ARTS_METHOD_ERROR_CATCH                            \
  catch (std::logic_error & e) {                           \
    throw arts::catch_errors(e, src_location{}.getfunc()); \
  }                                                        \
  catch (std::exception & e) {                             \
    throw arts::catch_errors(e, src_location{}.getfunc()); \
  }

#ifndef NDEBUG
#include <exception>

#else

#endif /* NDEBUG */

#endif /* debug_h */

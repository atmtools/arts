/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#include <format_tags.h>

#include <source_location>

struct src_location {
  std::source_location loc;
  src_location(std::source_location = std::source_location::current());
  std::string get();
  std::string getfunc();
};

namespace arts {
std::runtime_error catch_errors(std::logic_error& e,
                                const std::string_view context);
std::runtime_error catch_errors(std::exception& e,
                                const std::string_view context);

std::runtime_error user_error(const std::string_view msg,
                              const std::string_view context);
std::runtime_error user_error(const std::string_view msg,
                              const std::string_view condition,
                              const std::string_view context);
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

/*! An error has occured, will throw the error */
#define ARTS_USER_ERROR(fmt, ...)                                    \
  {                                                                  \
    throw arts::user_error(std::string_view{__VA_OPT__(std::format)( \
                               fmt __VA_OPT__(, ) __VA_ARGS__)},     \
                           src_location{}.get());                    \
  }

/*! Condition should be false to pass external check */
#define ARTS_USER_ERROR_IF(condition, fmt, ...)                        \
  {                                                                    \
    if (condition)                                                     \
      throw arts::user_error(std::string_view{__VA_OPT__(std::format)( \
                                 fmt __VA_OPT__(, ) __VA_ARGS__)},     \
                             #condition,                               \
                             src_location{}.get());                    \
  }

#endif /* debug_h */

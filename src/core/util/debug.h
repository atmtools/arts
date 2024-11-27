/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#include <cassert>
#include <format>
#include <version>

// These three overloads exist to allow for zero arguments to be passed to std::format
template <typename... Args>
std::string artsformat(std::format_string<Args...> fmt, Args&&... args)
  requires(sizeof...(Args) > 0)
{
  return std::format(fmt, std::forward<Args>(args)...);
}
std::string artsformat(std::format_string<> fmt);
std::string artsformat();

#if __cpp_lib_source_location >= 201907L
#include <source_location>

#define CURRENT_SOURCE_LOCATION                                          \
  std::format(                                                           \
      "Filename:      {}\n"                                              \
      "Function Name: {}\n"                                              \
      "Line Number:   {}\n"                                              \
      "Column Number: {}\n",                                             \
      std::string_view(std::source_location::current().file_name()),     \
      std::string_view(std::source_location::current().function_name()), \
      std::source_location::current().line(),                            \
      std::source_location::current().column())

#else

#define CURRENT_SOURCE_LOCATION \
  std::format("File:Line:     " __FILE__ ":{}", __LINE__)

#endif

#if __cpp_lib_source_location >= 201907L
#define CURRENT_SOURCE_FUNCTION \
  std::string(std::source_location::current().function_name())
#else
#define CURRENT_SOURCE_FUNCTION std::format(__FILE__ ":{}", __LINE__)
#endif

#define ARTS_METHOD_ERROR_CATCH                                                \
  catch (std::logic_error & e) {                                               \
    throw std::runtime_error(                                                  \
        std::format("Assertion error caught in:\n{}\n\n{}",                    \
                    CURRENT_SOURCE_FUNCTION,                                   \
                    std::string_view(e.what())));                              \
  }                                                                            \
  catch (std::exception & e) {                                                 \
    throw std::runtime_error(std::format(                                      \
        "{}\n{}", CURRENT_SOURCE_FUNCTION, std::string_view(e.what())));       \
  }

#ifndef NDEBUG
#include <exception>

#ifdef ARTS_ASSERT_USE_C

#define ARTS_ASSERT(condition, ...) \
  {                                 \
    assert(condition);              \
  }

#else

/*! Condition should be true to pass internal check */
#define ARTS_ASSERT(condition, ...)                                     \
  {                                                                     \
    if (not(condition)) {                                               \
      throw std::logic_error(                                           \
          std::format("Failed Internal Assert: {}\n"                    \
                      "This is a bug!  Bug is found at:\n{}\n"          \
                      "\nPlease contact ARTS developers so we can fix " \
                      "our error(s) via:\n\t"                           \
                      "github.com/atmtools/arts\n{}",                   \
                      std::string_view(#condition),                     \
                      CURRENT_SOURCE_LOCATION,                          \
                      artsformat(__VA_OPT__(__VA_ARGS__))));            \
    }                                                                   \
  }

#endif

#else

/*! Condition should be true to pass internal check, lets hope it is! */
#define ARTS_ASSERT(condition, ...) \
  {                                 \
  }

#endif /* NDEBUG */

/*! An error has occured, will throw the error */
#define ARTS_USER_ERROR(...)                                         \
  {                                                                  \
    static_assert(false __VA_OPT__(or true),                         \
                  "Must have an error message in user-"              \
                  "facing code in " __FILE__);                       \
    if constexpr (false __VA_OPT__(or true))                         \
      throw std::runtime_error(std::format(                          \
          "User error found at:\n"                                   \
          "\n"                                                       \
          "{}\n"                                                     \
          "\n"                                                       \
          "Please follow these instructions to correct the error:\n" \
          "\n"                                                       \
          "{}\n",                                                    \
          CURRENT_SOURCE_LOCATION,                                   \
          std::format(__VA_OPT__(__VA_ARGS__))));                    \
  }

/*! Condition should be false to pass external check */
#define ARTS_USER_ERROR_IF(condition, ...)                             \
  {                                                                    \
    if (condition) {                                                   \
      static_assert(false __VA_OPT__(or true),                         \
                    "Must have an error message in user-"              \
                    "facing code in " __FILE__);                       \
      if constexpr (false __VA_OPT__(or true))                         \
        throw std::runtime_error(std::format(                          \
            "Fail check: {}\n"                                         \
            "\n"                                                       \
            "User error found at:\n"                                   \
            "\n"                                                       \
            "{}\n"                                                     \
            "\n"                                                       \
            "Please follow these instructions to correct the error:\n" \
            "\n"                                                       \
            "{}\n",                                                    \
            std::string_view(#condition),                              \
            CURRENT_SOURCE_LOCATION,                                   \
            std::format(__VA_OPT__(__VA_ARGS__))));                    \
    }                                                                  \
  }

#endif /* debug_h */

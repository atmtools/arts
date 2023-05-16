/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#include <sstream>
#include <string>
#include <version>

/*! Take all arguments and turn to string by their operator<<() */
template <typename ... Args>
std::string var_string(Args&& ... args) {
  if constexpr (sizeof...(Args) not_eq 0) {
    std::ostringstream os;
    ((os << std::forward<Args>(args)), ...);
    return os.str();
  } else {
    return "";
  }
}

#if __cpp_lib_source_location	>= 201907L
#include <iomanip>
#include <source_location>

#define CURRENT_SOURCE_LOCATION                                                \
  var_string("Filename:      ",                                                \
             std::quoted(std::source_location::current().file_name()), '\n',   \
             "Function Name: ",                                                \
             std::quoted(std::source_location::current().function_name()),     \
             '\n', "Line Number:   ", std::source_location::current().line(),  \
             '\n',                                                             \
             "Column Number: ", std::source_location::current().column())

#else

#define CURRENT_SOURCE_LOCATION var_string("\t", __FILE__, ":", __LINE__)

#endif

#ifndef NDEBUG
#include <exception>
#include <iostream>

// Use this macro around function parameter names and variable definitions
// which are only used in assertions
#define DEBUG_ONLY(...) __VA_ARGS__

// Use this macro to output a counter value everytime a
// certain place is reached
#define DEBUG_COUNTER(n)                                    \
  {                                                         \
    static Index n = 0;                                     \
    std::cerr << "DBG: " << #n << ": " << ++n << std::endl; \
  }

// Print value of expression for debugging
#define DEBUG_PRINT(e) \
  { std::cerr << "DBG: " << (e) << std::endl; }

// Print expression and value for debugging
#define DEBUG_VAR(e) \
  { std::cerr << "DBG: " << #e << ": " << (e) << std::endl; }

// Print expression and value with the given fp precision for debugging
#define DEBUG_VAR_FLT(p, e)                                           \
  {                                                                   \
    std::streamsize old_p = std::cerr.precision();                    \
    std::cerr << "DBG: " << #e << ": " << std::setprecision(p) << (e) \
              << std::endl                                            \
              << std::setprecision(old_p);                            \
  }

/*! Turn off noexcept */
#define ARTS_NOEXCEPT noexcept(false)

/*! Condition should be true to pass internal check */
#define ARTS_ASSERT(condition, ...)                                            \
  {                                                                            \
    if (not(condition)) {                                                      \
      throw std::runtime_error(                                                \
          var_string("Failed Internal Assert: " #condition "\n"                \
                     "This is a bug!  Bug is found at:\n",                     \
                     CURRENT_SOURCE_LOCATION,                                  \
                     "\nPlease contact"                                        \
                     " ARTS developers so we can fix "                         \
                     "our error(s) via:\n\t"                                   \
                     "github.com/atmtools/arts\n" __VA_OPT__(, __VA_ARGS__))); \
    }                                                                          \
  }

#else

#define DEBUG_ONLY(...)

#define DEBUG_COUNTER(n)

#define DEBUG_PRINT(e)

#define DEBUG_VAR(e)

#define DEBUG_VAR_FLT(p, e)

/*! Turn on noexcept explicitly */
#define ARTS_NOEXCEPT noexcept(true)

/*! Condition should be true to pass internal check, lets hope it is! */
#define ARTS_ASSERT(condition, ...) {}

#endif /* NDEBUG */

#ifdef NO_ARTS_USER_ERRORS

/*! Turn on noexcept for user-facing functions */
#define ARTS_USER_NOEXCEPT noexcept(true)

/*! Condition should be false to pass external check, lets hope it is!  */
#define ARTS_USER_ERROR_IF(condition, ...) {}

/*! An error has occured, should throw the error, but will instead ignore */
#define ARTS_USER_ERROR(...) {}

#else  // NO_ARTS_USER_ERRORS == 0

/*! Turn off noexcept for user-facing functions */
#define ARTS_USER_NOEXCEPT noexcept(false)

/*! Condition should be false to pass external check */
#define ARTS_USER_ERROR_IF(condition, ...)                                     \
  {                                                                            \
    static_assert(false __VA_OPT__(or true),                                   \
                  "Must have an error message in user-"                        \
                  "facing code in " __FILE__);                                 \
    if (condition) {                                                           \
      throw std::runtime_error(                                                \
          var_string("User Error: " #condition "\nError is found at:\n",       \
                     CURRENT_SOURCE_LOCATION,                                  \
                     "\nPlease follow "                                        \
                     "these instructions to correct your "                     \
                     "error:\n" __VA_OPT__(, __VA_ARGS__)));                   \
    }                                                                          \
  }

/*! An error has occured, will throw the error */
#define ARTS_USER_ERROR(...)                                                   \
  {                                                                            \
    static_assert(false __VA_OPT__(or true),                                   \
                  "Must have an error message in user-"                        \
                  "facing code in " __FILE__);                                 \
    throw std::runtime_error(                                                  \
        var_string("User Error:\n"                                             \
                   "Error is found at:\n",                                     \
                   CURRENT_SOURCE_LOCATION,                                    \
                   "\nPlease follow "                                          \
                   "these instructions to correct your "                       \
                   "error:\n" __VA_OPT__(, __VA_ARGS__)));                     \
  }

#endif  // NO_ARTS_USER_ERRORS

#endif /* debug_h */

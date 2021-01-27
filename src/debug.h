/* Copyright (C) 2013
   Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#include <sstream>
#include <string>
#include <tuple>

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
#define ARTS_NOEXCEPT

/*! Condition should be true to pass internal check */
#define ARTS_ASSERT(condition, ...) {     \
  if (not (condition)) {                  \
    throw std::runtime_error(             \
    var_string("Failed Internal Assert: " \
        #condition "\n" "This is a bug! " \
        "Bug is found at:\n\t" __FILE__   \
        ":", __LINE__, "\nPlease contact" \
        " ARTS developers so we can fix " \
        "our error(s) via:\n\t"           \
        "github.com/atmtools/arts\n") +   \
    var_string(__VA_ARGS__)               \
    );                                    \
  } }

#else

#define DEBUG_ONLY(...)

#define DEBUG_COUNTER(n)

#define DEBUG_PRINT(e)

#define DEBUG_VAR(e)

#define DEBUG_VAR_FLT(p, e)

/*! Turn on noexcept explicitly */
#define ARTS_NOEXCEPT noexcept

/*! Condition should be true to pass internal check, lets hope it is! */
#define ARTS_ASSERT(condition, ...) {}

#endif /* NDEBUG */

#if NO_ARTS_USER_ERRORS == 1

/*! Turn on noexcept for user-facing functions */
#define ARTS_USER_NOEXCEPT noexcept

/*! Condition should be false to pass external check, lets hope it is!  */
#define ARTS_USER_ERROR_IF(condition, ...) {}

#else  // NO_ARTS_USER_ERRORS == 0

/*! Turn off noexcept for user-facing functions */
#define ARTS_USER_NOEXCEPT

/*! Condition should be false to pass external check */
#define ARTS_USER_ERROR_IF(condition, ...) {  \
  static_assert(std::tuple_size<decltype(     \
    std::make_tuple(__VA_ARGS__))>::value,    \
    "Must have an error message in user-"     \
    "facing code in " __FILE__);              \
  if (condition) {                            \
    throw std::runtime_error(                 \
      var_string("User Error: " #condition    \
        "\nError is found at:\n\t" __FILE__   \
        ":", __LINE__, "\nPlease follow "     \
        "these instructions to correct your " \
        "error:\n") + var_string(__VA_ARGS__) \
    );                                        \
  } }

#endif  // NO_ARTS_USER_ERRORS

#endif /* debug_h */

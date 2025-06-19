#ifndef depr_h
#define depr_h

#include <ios>
#include <iostream>
#include <streambuf>

#include "debug.h"

#ifndef NO_DEPR
#define DEPRECATED_FUNCTION(FUNCTION_NAME, DATE_STRING_OR_SILLY_REASON, ...) \
  static_assert(                                                             \
      false __VA_OPT__(or true),                                             \
      "Must give a deprecation reason in user facing code in " __FILE__);    \
  {                                                                          \
    static bool i_have_not_informed_you = true;                              \
    if (i_have_not_informed_you) {                                           \
      std::cerr << "# DEPRECATED FUNCTION WARNING\n"                         \
                << '\n'                                                      \
                << '\n'                                                      \
                << FUNCTION_NAME << " is deprecated since "                  \
                << DATE_STRING_OR_SILLY_REASON << '\n'                       \
                << '\n'                                                      \
                << std::format(__VA_ARGS__) << '\n'                          \
                << '\n'                                                      \
                << "# /DEPRECATED FUNCTION WARNING\n";                       \
      i_have_not_informed_you = false;                                       \
    }                                                                        \
  }
#endif

#endif  // depr_h

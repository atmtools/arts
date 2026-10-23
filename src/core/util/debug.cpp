#include "debug.h"

#include <mystring.h>

#include <exception>
#include <format>
#include <string>
#include <string_view>

std::string src_location::get() const {
  return std::format(
      "Filename:      {}\n"
      "Function Name: {}\n"
      "Line Number:   {}\n"
      "Column Number: {}\n",
      std::string_view(file_name_),
      std::string_view(function_name_),
      line_,
      column_);
}

std::string src_location::getfunc() const { return function_name_; }

namespace arts {
std::runtime_error catch_errors(std::logic_error& e, const std::string_view context) {
  return std::runtime_error(std::format("Assertion error caught in:\n{}\n\n{}", context, indent(e.what())));
}

std::runtime_error catch_errors(std::exception& e, const std::string_view context) {
  return std::runtime_error(std::format("{}\n{}", context, indent(e.what())));
}

std::runtime_error user_error(const std::string_view msg, const src_location& location) {
#ifdef NDEBUG
  return std::runtime_error(std::format(R"(Function Name: {1}

{0})",
                                        msg,
                                        location.getfunc()));
#else
  return std::runtime_error(std::format(
      R"(User error found in:                                   
                                                     
{1}
                                                     
{0})",
      msg,
      location.get()));
#endif
}

std::runtime_error user_error(const std::string_view                  msg,
                              [[maybe_unused]] const std::string_view condition,
                              const src_location&                     location) {
#ifdef NDEBUG
  return user_error(msg, location);
#else
  return std::runtime_error(std::format(
      R"(Fail check: {1}

User error found in:                                   
                                                     
{2}   
                                                    
Please follow these instructions to correct the error:
                                                     
{0})",
      msg,
      condition,
      location.get()));
#endif
}
}  // namespace arts
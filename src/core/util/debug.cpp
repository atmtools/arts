#include "debug.h"

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
  return std::runtime_error(std::format("Assertion error caught in:\n{}\n\n{}", context, std::string_view(e.what())));
}

std::runtime_error catch_errors(std::exception& e, const std::string_view context) {
  return std::runtime_error(std::format("{}\n{}", context, std::string_view(e.what())));
}

std::runtime_error user_error(const std::string_view msg, const std::string_view context) {
  return std::runtime_error(std::format(
      R"(User error found in:                                   
                                                     
{1}   
                                                    
Please follow these instructions to correct the error:
                                                     
{0})",
      msg,
      context));
}

std::runtime_error user_error(const std::string_view msg,
                              const std::string_view condition,
                              const std::string_view context) {
  return std::runtime_error(std::format(
      R"(Fail check: {1}

User error found in:                                   
                                                     
{2}   
                                                    
Please follow these instructions to correct the error:
                                                     
{0})",
      msg,
      condition,
      context));
}
}  // namespace arts
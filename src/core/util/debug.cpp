#include "debug.h"

#include <exception>
#include <format>
#include <string>
#include <string_view>

#if __cpp_lib_source_location >= 201907L
src_location::src_location(std::source_location loc_) : loc(loc_) {}
#else
src_location::src_location(std::string_view f_, const int l_) : f(f_), l(l_) {}
#endif

std::string src_location::get() {
  return
#if __cpp_lib_source_location >= 201907L
      std::format(
          "Filename:      {}\n"
          "Function Name: {}\n"
          "Line Number:   {}\n"
          "Column Number: {}\n",
          std::string_view(loc.file_name()),
          std::string_view(loc.function_name()),
          loc.line(),
          loc.column());
#else
      std::format("File:Line:     {}:{}", f, l);
#endif
}

std::string src_location::getfunc() {
  return
#if __cpp_lib_source_location >= 201907L
      std::source_location::current().function_name();
#else
      get();
#endif
}

namespace arts {
std::runtime_error catch_errors(std::logic_error& e,
                                const std::string_view context) {
  return std::runtime_error(std::format("Assertion error caught in:\n{}\n\n{}",
                                        context,
                                        std::string_view(e.what())));
}

std::runtime_error catch_errors(std::exception& e,
                                const std::string_view context) {
  return std::runtime_error(
      std::format("{}\n{}", context, std::string_view(e.what())));
}

std::runtime_error user_error(const std::string_view msg,
                              const std::string_view context) {
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
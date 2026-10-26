#include "debug.h"

#include <mystring.h>

#include <atomic>
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

namespace {
std::atomic<bool> show_cpp_error_location_ {
#ifdef NDEBUG
  false
#else
  true
#endif
};
}  // namespace

bool show_cpp_error_location() { return show_cpp_error_location_.load(std::memory_order_relaxed); }

void show_cpp_error_location(bool show) { show_cpp_error_location_.store(show, std::memory_order_relaxed); }

void user_error(const std::string_view msg, const src_location& location) {
  if (not show_cpp_error_location()) throw std::runtime_error(std::string{msg});

  throw std::runtime_error(std::format("ARTS C++ error location:\n{}\n{}", indent(location.get()), msg));
}
}  // namespace arts
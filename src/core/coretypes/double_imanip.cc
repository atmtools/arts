#include "double_imanip.h"

#include <debug.h>
#include <fast_float/fast_float.h>

#include <istream>
#include <system_error>

const double_imanip& double_imanip::operator>>(double& x) const {
  std::istream& is = *in;
  std::string   buf;

  // Read to the buffer
  is >> buf;
  if (is.fail()) arts::user_error("Cannot read from stream");

  // Actual conversion
  const auto res = fast_float::from_chars(buf.c_str(), buf.c_str() + buf.size(), x);

  // Error (only std::errc::invalid_argument possible)
  if (res.ec == std::errc::invalid_argument)
    arts::user_error(R"--(The argument:

`{}`

is not convertible to a valid double.  At the very least it
cannot be converted to one using the standard string-to-double
routine
)--",
                     buf);

  if (!is.eof() && is.tellg() != -1) {
    is.seekg(std::distance(buf.c_str(), res.ptr) - buf.size(), std::ios_base::cur);
  } else {
    // Use alternative method for non-seekable files
    // This is known to break in very rare cases on macOS with Clang
    std::size_t n = std::distance(buf.c_str(), res.ptr);
    while (n++ < buf.size()) is.unget();
  }

  return *this;
}

std::istream& double_imanip::operator>>(const double_imanip&) const { return *in; }

const double_imanip& operator>>(std::istream& in, const double_imanip& dm) {
  dm.in = &in;
  return dm;
}
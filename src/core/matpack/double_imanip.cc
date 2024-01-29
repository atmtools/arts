#include "double_imanip.h"

#include <fast_float/fast_float.h>

#include <ios>
#include <system_error>

#include "debug.h"

const double_imanip& double_imanip::operator>>(double& x) const {
  std::istream& is = *in;
  std::string buf;

  // Read to the buffer
  is >> buf;
  ARTS_USER_ERROR_IF(is.fail(), "Cannot read from stream")

  // Actual conversion
  const auto res =
      fast_float::from_chars(buf.c_str(), buf.c_str() + buf.size(), x);

  // Error (only std::errc::invalid_argument possible)
  ARTS_USER_ERROR_IF(res.ec == std::errc::invalid_argument,
                     "The argument: \n\n'",
                     buf,
                     R"--('

is not convertible to a valid double.  At the very least it
cannot be converted to one using the standard string-to-double
routine
)--")

  if (!is.eof() && is.tellg() != -1) {
    is.seekg(std::distance(buf.c_str(), res.ptr) - buf.size(),
             std::ios_base::cur);
  } else {
    // Use alternative method for non-seekable files
    // This is known to break in very rare cases on macOS with Clang
    std::size_t n = std::distance(buf.c_str(), res.ptr);
    while (n++ < buf.size()) is.unget();
  }

  return *this;
}
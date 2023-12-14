#pragma once

#include <fast_float/fast_float.h>

#include <charconv>
#include <sstream>
#include <string>
#include <string_view>

#include "array.h"
#include "configtypes.h"

/** The String type for ARTS. Implementation. */
using String = std::string;

/** An array of Strings. */
using ArrayOfString = Array<String>;

/** An array of Strings. */
using ArrayOfArrayOfString = Array<Array<String>>;

/** Name string_view as we named string */
using StringView = std::string_view;

namespace std {
std::ostream& operator<<(std::ostream& os, const ArrayOfString& x);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfString& x);
}  // namespace std

/** Extract something from the beginning of a string. This is just a small helper
 function to safe some typing.

 \retval x    What was extracted from the beginning of the line.
 \retval line What was extracted is also cut away from line.
 \param n     The width of the stuff to extract.

 \author Stefan Buehler */
template <class T>
void extract(T& x, String& line, Size n) {
  // Initialize output to zero! This is important, because otherwise
  // the output variable could `remember' old values.
  x = T(0);

  const Size N = n;
  Size i = 0;
  while (i < N and i < line.size() and isspace(line[i])) ++i;
  while (n > i and (n - 1) < line.size() and isspace(line[n - 1])) --n;

  if constexpr (std::is_same_v<double, T> or std::is_same_v<float, T>) {
    fast_float::from_chars(line.data() + i, line.data() + n, x);
  } else if constexpr (std::is_same_v<long long, T> or
                       std::is_same_v<long, T> or std::is_same_v<int, T>) {
    std::from_chars(line.data() + i, line.data() + n, x);
  } else {
    // This will contain the short subString with the item to extract.
    // Make it a String stream, for easy parsing,
    // extracting subString of width n from line:
    std::istringstream item(line.substr(i, n));

    // Convert with the aid of String stream item:
    item >> x;
  }

  // Shorten line by n:
  line.erase(0, N);
}

void tolower(String& x);
String tolower(const String& x);

void toupper(String& x);
String toupper(const String& x);

void split(ArrayOfString& aos, const String& x, const String& delim);
ArrayOfString split(const String& x, const String& delim);

void trim(String& x);
String trim(const String& x);

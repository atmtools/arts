#pragma once

#include <fast_float/fast_float.h>

#include <algorithm>
#include <charconv>
#include <string>
#include <string_view>

#include "array.h"
#include "nonstd.h"
#include <matpack_concepts.h>

/** The String type for ARTS. Implementation. */
using String = std::string;

/** An array of Strings. */
using ArrayOfString = Array<String>;

/** An array of Strings. */
using ArrayOfArrayOfString = Array<Array<String>>;

namespace std {
inline std::ostream& operator<<(std::ostream& os, const ArrayOfString& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}

inline std::ostream& operator<<(std::ostream& os,
                                const ArrayOfArrayOfString& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}
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
  while (n > i and (n-1) < line.size() and isspace(line[n-1])) --n;
  
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

inline
void tolower(String& x) {
  std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c){ return ::tolower(c); });
}

[[nodiscard]]
inline
String tolower(const String& x) {
  String out=x;
  tolower(out);
  return out;
}

inline
void toupper(String& x) {
  std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c){ return ::toupper(c); });
}

[[nodiscard]]
inline
String toupper(const String& x) {
  String out=x;
  toupper(out);
  return out;
}

inline
void split(ArrayOfString& aos, const String& x, const String& delim) {
    Size pos, oldpos;
    pos = oldpos = 0;
    aos.resize(0);

    while (oldpos < x.size() &&
           (pos = x.find(delim, oldpos)) != x.npos) {
      if (pos && pos - oldpos)
        aos.push_back(x.substr(oldpos, pos - oldpos));
      oldpos = pos + delim.size();
    }

    if (oldpos < x.size()) aos.push_back(x.substr(oldpos));
}

[[nodiscard]]
inline
ArrayOfString split(const String& x, const String& delim) {
  ArrayOfString out;
  split(out, x, delim);
  return out;
}

inline
void trim(String& x) {
  while (nonstd::isspace(x.front())) x.erase(x.begin());
  while (nonstd::isspace(x.back())) x.pop_back();
}

[[nodiscard]]
inline
String trim(const String& x) {
  String out = x;
  trim(out);
  return out;
}

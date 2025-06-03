#pragma once

#include "string_extract.h"

#include <fast_float/fast_float.h>

#include <charconv>
#include <sstream>
#include <string>

/** Extract something from the beginning of a string. This is just a small helper
 function to safe some typing.

 \retval x    What was extracted from the beginning of the line.
 \retval line What was extracted is also cut away from line.
 \param n     The width of the stuff to extract.

 \author Stefan Buehler */
template <class T>
void extract_tmpl(T& x, std::string& line, std::size_t n) {
  // Initialize output to zero! This is important, because otherwise
  // the output variable could `remember' old values.
  x = T(0);

  const std::size_t N = n;
  std::size_t i       = 0;
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

void extract(float& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(double& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(char& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(int& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(long int& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(long long int& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(unsigned char& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(unsigned int& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(unsigned long int& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}
void extract(unsigned long long int& x, std::string& line, std::size_t n){extract_tmpl(x, line, n);}

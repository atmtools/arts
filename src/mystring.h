/* Copyright (C) 2001-2012 Stefan Buehler <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/** 
    \file   mystring.h

   This file contains the definition of String, the ARTS string class.

   \author Stefan Buehler
   \date   2001-09-14
*/

#ifndef mystring_h
#define mystring_h

#include <algorithm>
#include <climits>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include "array.h"
#include "matpack.h"
#include "double_imanip.h"

/**
   The implementation for String, the ARTS string class. 

   This adds some additional functionality to the standard stl string
   class, notably:

   a) Range checking by assert

   b) nelem() member function, return the size of the String of type
   Index. 

   The type string is just a typedef for
   basic_string<char>. Therefore, to make everything work
   correctly, we have to derive our own class from basic_string,
   not from string directly.
*/
template <class charT>
class my_basic_string : public std::basic_string<charT> {
 public:
  // Constructors:
  my_basic_string() = default;
  explicit my_basic_string(Index n, char c = ' ');
  my_basic_string(const std::basic_string<charT>& A,
                  Index pos = 0,
                  Index numpos = my_basic_string<charT>::npos);
  my_basic_string(const char A[]);
  my_basic_string(const std::string_view& sv);

  // Insert string before all occurrences of the substring.
  void insert_substr(const my_basic_string<charT>& searchstr,
                     const my_basic_string<charT>& insstr);

  // Split string
  void split(Array<my_basic_string<charT> >& aos,
             const my_basic_string<charT>& delim) const;

  /** Convert to upper case */
  void toupper() {
    std::transform(this->begin(), this->end(), this->begin(), ::toupper);
  }

  my_basic_string toupper() const {
    my_basic_string s = *this;
    s.toupper();
    return s;
  }

  /** Convert to lower case */
  void tolower() {
    std::transform(this->begin(), this->end(), this->begin(), ::tolower);
  }

  my_basic_string tolower() const {
    my_basic_string s = *this;
    s.tolower();
    return s;
  }

  /** Trim leading and trailing whitespace */
  void trim();

  // Number of elements:
  Index nelem() const;

  // Index operators:
  char operator[](Index n) const;
  char& operator[](Index n);

  /** Define npos: */
  static const Index npos = static_cast<Index>(std::basic_string<charT>::npos);

  typedef Index size_type;
};

// Member functions for my_basic_string:

// Constructors:

/** Constructor setting size. You may give as a second argument a
    character with which to fill the new string. Per default this is
    zero. 
    
    \param n Number of characters
    \param c Optional fill character
*/
template <class charT>
inline my_basic_string<charT>::my_basic_string(Index n, char c)
    : std::basic_string<charT>(n, c) { /* Nothing to do here. */
}

/** Construnctor from a basic_string. This is important for handling
    of expressions like this to work correctly:

    String a = b+'.'+c 

    As for basic_string, this constructor can also be used to
    initialize the new string from a subrange of the original string. 

    \param A The original string
    \param pos Start position (0 means from the beginning)
    \param numpos How many characters to copy
*/
template <class charT>
inline my_basic_string<charT>::my_basic_string(
    const std::basic_string<charT>& A, Index pos, Index numpos) {
  // Range checks:
  ARTS_ASSERT(0 <= pos);  // Start index must be 0 or greater 0.

  if (!A.size()) return;

  //   cout << "A = " << A << "\n";
  //   cout << "pos = " << pos << "\n";
  //   cout << "size = " << A.size() << "\n";

  ARTS_ASSERT(static_cast<typename std::basic_string<charT>::size_type>(pos) <
         A.size());
  // At most the last element of the original string.

  ARTS_ASSERT(numpos == my_basic_string<charT>::npos ||
         ((numpos >= 0) &&
          (static_cast<typename std::basic_string<charT>::size_type>(numpos) <=
           (A.size() -
            pos))));  // Number of characters to copy must be at the most the
                      // number left. -1 means all remaining characters.

  // The assertions look complicated, because we have to cast pos and
  // npos to the unsigned size type of basic string to avoid warning
  // messages from the compiler. Both casts are save, because previous
  // assertions check that pos and npos are positive. (The allowed
  // case npos -1 (=my_basic_string<charT>::npos) is also handled
  // correctly.)

  std::basic_string<charT>::operator=(std::basic_string<charT>(A, pos, numpos));
}

/** Constructor from a C-style char array. */
template <class charT>
inline my_basic_string<charT>::my_basic_string(const char A[])
    : std::basic_string<charT>(A) { /* Nothing to do here. */
}

/** Constructor from a std::string_view. */
template <class charT>
inline my_basic_string<charT>::my_basic_string(const std::string_view& sv)
    : std::basic_string<charT>(std::string(sv)) { /* Nothing to do here. */
}

/** Insert string before all occurrences of the substring.
 
 \param[in] searchstr  String to search for.
 \param[in] insstr     String to insert.
*/
template <class charT>
inline void my_basic_string<charT>::insert_substr(
    const my_basic_string<charT>& searchstr,
    const my_basic_string<charT>& insstr) {
  size_t searchstr_size = searchstr.size();
  size_t insstr_size = insstr.size();
  size_t start_pos = 0;

  while (start_pos != std::string::npos) {
    start_pos = this->find(searchstr, start_pos);
    if (start_pos && start_pos != std::string::npos) {
      this->insert(start_pos, insstr);
      start_pos += searchstr_size + insstr_size;
    }
  }
}

/** Split string into substrings.
 
 \param[out] aos    ArrayOfString containing the returned substrings.
 \param[in]  delim  Delimiter string.
 */
template <class charT>
inline void my_basic_string<charT>::split(
    Array<my_basic_string<charT> >& aos,
    const my_basic_string<charT>& delim) const {
  size_t pos, oldpos;
  pos = oldpos = 0;
  aos.resize(0);

  while (oldpos < (size_t)this->nelem() &&
         (pos = this->find(delim, oldpos)) !=
             (size_t)my_basic_string<charT>::npos) {
    if (pos && pos - oldpos) aos.push_back(this->substr(oldpos, pos - oldpos));
    oldpos = pos + delim.nelem();
  }

  if (oldpos < (size_t)this->nelem()) aos.push_back(this->substr(oldpos));
}

/** Trim leading and trailing whitespace */
template <class charT>
void my_basic_string<charT>::trim() {
  // Create ref to self for readability
  my_basic_string& this_string = *this;

  // Remove leading whitespace
  while (0 != this_string.nelem() &&
         (' ' == this_string[0] || '\t' == this_string[0] ||
          '\n' == this_string[0] || '\r' == this_string[0]))
    this_string.erase(0, 1);

  // Remove trailing whitespace
  while (0 != this_string.nelem() &&
         (' ' == this_string[this_string.nelem() - 1] ||
          '\t' == this_string[this_string.nelem() - 1] ||
          '\n' == this_string[this_string.nelem() - 1] ||
          '\r' == this_string[this_string.nelem() - 1]))
    this_string.erase(this_string.nelem() - 1);
}

/** Number of elements. */
template <class charT>
inline Index my_basic_string<charT>::nelem() const {
  size_t s = this->size();
  ARTS_ASSERT(s < LONG_MAX);
  return static_cast<long>(s);
}

/**
  Constant index operator. We redifine this here so that we can have
  range checking by assert.
    
  \param[in] n Index  
*/
template <class charT>
inline char my_basic_string<charT>::operator[](Index n) const {
  ARTS_ASSERT(0 <= n);
  ARTS_ASSERT(n < nelem());
  return std::basic_string<charT>::operator[](n);
}

/**
  Non-constant index operator. We redifine this here so that we can
  have range checking by assert.
    
  \param[in] n Index
*/
template <class charT>
inline char& my_basic_string<charT>::operator[](Index n) {
  ARTS_ASSERT(0 <= n);
  ARTS_ASSERT(n < nelem());
  return std::basic_string<charT>::operator[](n);
}

/** The String type for ARTS. Implementation see documentation of
    class my_basic_string. */
typedef my_basic_string<char> String;

/** An array of Strings. */
typedef Array<String> ArrayOfString;

/** An array of Strings. */
typedef Array<Array<String> > ArrayOfArrayOfString;

/** Extract something from the beginning of a string. This is just a small helper
 function to safe some typing.

 \retval x    What was extracted from the beginning of the line.
 \retval line What was extracted is also cut away from line.
 \param n     The width of the stuff to extract.

 \author Stefan Buehler */
template <class T>
void extract(T& x, String& line, Index n) {
  // Initialize output to zero! This is important, because otherwise
  // the output variable could `remember' old values.
  x = T(0);

  if constexpr (std::is_same_v<double, T> or std::is_same_v<float, T>) {
    fast_float::from_chars(line.data(), line.data() + n, x);
  } else {
    // This will contain the short subString with the item to extract.
    // Make it a String stream, for easy parsing,
    // extracting subString of width n from line:
    std::istringstream item(line.substr(0, n));

    //  cout << "line = '" << line << "'\n";
    //   cout << "line.substr(0,n) = " << line.substr(0,n) << endl;
    //   cout << "item = " << item.str() << endl;
    //  cout << "line = " << line << endl;

    // Convert with the aid of String stream item:
    item >> x;
  }

  // Shorten line by n:
  line.erase(0, n);
}

#endif  // mystring_h

// <!!-*-C++-*- file: iotraits.h --->
// <!!-------------------------------------------------------------------------->
// <!! Copyright (C) 1996 Dietmar Kuehl >
// <!!   Universitaet Konstanz, Lehrstuhl fuer praktische Informatik I >
// <!!>
// <!! Permission to use, copy, modify, distribute and sell this >
// <!! software for any purpose is hereby granted without fee, provided >
// <!! that the above copyright notice appears in all copies and that >
// <!! both that copyright notice and this permission notice appear in >
// <!! supporting documentation. Dietmar Kuehl makes no representations >
// <!! about the suitability of this software for any purpose. It is >
// <!! provided "as is" without express or implied warranty. >
// <!!-------------------------------------------------------------------------->

// Author: Dietmar Kühl dietmar.kuehl@uni-konstanz.de www.informatik.uni-konstanz.de/~kuehl
// Title:  The definition of the iotraits classes

//------------------------------------------------------------------------------

#if !defined(IOTRAITS_H)
#define IOTRAITS_H

#include <cstring>
#include <cctype>
#include <streambuf.h>

#if defined(__GNUC__) && __GNUC_MINOR__ <= 7
#  define PORT_TEMPLATE_SPECIALIZATION
#else
#  define PORT_TEMPLATE_SPECIALIZATION template <>
#endif

//------------------------------------------------------------------------------

template <class charT> struct ctype {};

PORT_TEMPLATE_SPECIALIZATION
struct ctype<char>
{
  typedef char         char_type;
  typedef int          int_type;

  //--- fix needed:  bool      isspace(int_type c) const { return isspace(c); }
  char_type newline() const           { return '\n'; }
};

//------------------------------------------------------------------------------

template <class charT> struct ios_traits {};

//------------------------------------------------------------------------------

PORT_TEMPLATE_SPECIALIZATION
struct ios_traits<char>
{
  typedef char      char_type;
  typedef int       int_type;
  typedef streampos pos_type;
  typedef streamoff off_type;
  typedef int       state_type;

  static inline char_type eos()                { return char_type(); }
  static inline int_type  eof()                { return -1; }
  static inline int_type  not_eof(int_type c)  { return c == eof()? eof() + 1: c; }
  static inline char_type newline()            { return '\n'; }

  static inline bool eq(char_type c1, char_type c2) { return c1 == c2; }
  static inline bool eq_int_type(int_type i1, int_type i2)    { return i1 == i2; }
  static inline void assign(char& to, char from) { to = from; }

  static inline char_type to_char_type(int_type c) { return char_type(c); }
  static inline int_type  to_int_type(char_type c) { return int_type(c); }

  static inline char_type *copy(char_type *dst, char_type const *src, pos_type n)
    {
      memcpy(dst, src, n);
      return dst;
    }
  static inline size_t    length(char_type const *s) { return strlen(s); }
};

//------------------------------------------------------------------------------

#endif /* IOTRAITS_H */

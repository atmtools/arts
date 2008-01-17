/* Copyright (C) 2001-2007 Stefan Buehler <sbuehler@ltu.se>

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

#ifndef string_h
#define string_h

#include <cassert>
#include <climits>
#include <string>
#include "matpack.h"

// String stream library. This is included with the ARTS source code
// for now, because it is missing in gcc <= 2.95.2
#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include "sstream.h"
#endif

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
template<class charT>
class my_basic_string : public basic_string<charT>
{
public:
  // Constructors:
  my_basic_string();
  explicit my_basic_string(Index n, char c=' ');
  my_basic_string(const basic_string<charT>& A,
                  Index pos=0,
                  Index numpos=my_basic_string<charT>::npos);
  my_basic_string(const char A[]);

  // Assignment operators:
  my_basic_string& operator=(const my_basic_string<charT>& A);
  //  my_basic_string& operator=(const char A[]);

  // Number of elements:
  Index nelem() const;

  // Find functions:
//   Index find(char c);
//   Index find(const my_basic_string<charT>&  c);

  // Index operators:
  char operator[](Index n) const;
  char& operator[](Index n);

  /** Define npos: */
  static const Index npos = static_cast<Index>(basic_string<charT>::npos);

  typedef Index size_type;
};


// Member functions for my_basic_string:


// Constructors:

/** Default constructor. */
template<class charT>
inline my_basic_string<charT>::my_basic_string() : basic_string<charT>()  
{ /* Nothing to do here. */ }

/** Constructor setting size. You may give as a second argument a
    character with which to fill the new string. Per default this is
    zero. 
    
    \param n Number of characters
    \param c Optional fill character
*/
template<class charT>
inline my_basic_string<charT>::my_basic_string(Index n, char c) :
  basic_string<charT>(n,c) 
{ /* Nothing to do here. */ }

/** Construnctor from another my_basic_string. */
// template<class charT>
// inline my_basic_string<charT>::my_basic_string(const my_basic_string& A) : basic_string<charT>(A) 
// { /* Nothing to do here. */ };

/** Construnctor from a basic_string. This is important for handling
    of expressions like this to work correctly:

    String a = b+'.'+c 

    As for basic_string, this constructor can also be used to
    initialize the new string from a subrange of the original string. 

    \param A The original string
    \param pos Start position (0 means from the beginning)
    \param numpos How many characters to copy
*/
template<class charT>
inline my_basic_string<charT>::my_basic_string(const basic_string<charT>& A,
                                               Index pos,
                                               Index numpos)
{ 
  // Range checks:
  assert(0<=pos);               // Start index must be 0 or greater 0.

//   cout << "A = " << A << "\n";
//   cout << "pos = " << pos << "\n";
//   cout << "size = " << A.size() << "\n";

  assert(static_cast<typename basic_string<charT>::size_type>(pos)<A.size());      
  // At most the last element of the original string.

  assert( numpos==my_basic_string<charT>::npos ||
          ( (numpos >= 0) &&
            (static_cast<typename basic_string<charT>::size_type>(numpos)<=(A.size()-pos))
            )
          );  // Number of characters to copy must be at the most the
              // number left. -1 means all remaining characters. 

  // The assertions look complicated, because we have to cast pos and
  // npos to the unsigned size type of basic string to avoid warning
  // messages from the compiler. Both casts are save, because previous
  // assertions check that pos and npos are positive. (The allowed
  // case npos -1 (=my_basic_string<charT>::npos) is also handled
  // correctly.)

  basic_string<charT>::operator=(basic_string<charT>(A,pos,numpos));

}

/** Constructor from a C-style char array. */
template<class charT>
inline my_basic_string<charT>::my_basic_string(const char A[]) : basic_string<charT>(A) 
{ /* Nothing to do here. */ }


/** Assignment from another my_basic_string.

    The two partners do not have to have the same size. Size of the
    target string is adjusted automatically, just as string does
    it.
*/
template<class charT>
inline my_basic_string<charT>& my_basic_string<charT>::operator=(const my_basic_string<charT>& A)
{
  basic_string<charT>::operator=(A);
  return *this;
}

/** Number of elements. */
template<class charT>
inline Index my_basic_string<charT>::nelem() const
{ 
  size_t s = this->size();
  assert(s<LONG_MAX);
  return static_cast<long>(s);
}

// /** Find function for char.

//     \param c What character to find.
//     \return Position of c, or npos if not found.

//     Unfortunately, the basid_string.find() functions returns npos
//     when the character is not found. This is -1, but assigned to a
//     positive type! Gives a very high number. If you add 1 to this
//     number you get zero again, so it does indeed bahave like -1 in a
//     way. Yuck! With Index this does not work. */
// template<class charT>
// inline Index my_basic_string<charT>::find(char c)
// {
//   basic_string<charT>::size_type i = basic_string<charT>::find(c);
//   if ( i == basic_string<charT>::npos )
//     return npos;
//   else
//     return static_cast<Index>(i);
// }

// /** Find function for string.

//     \param c What string to find.
//     \return Position of c, or npos if not found.

//     Unfortunately, the basid_string.find() functions returns npos
//     when the character is not found. This is -1, but assigned to a
//     positive type! Gives a very high number. If you add 1 to this
//     number you get zero again, so it does indeed bahave like -1 in a
//     way. Yuck! With Index this does not work. */
// template<class charT>
// inline Index my_basic_string<charT>::find(const my_basic_string<charT>& c)
// {
//   basic_string<charT>::size_type i = basic_string<charT>::find(c);
//   if ( i == basic_string<charT>::npos )
//     return npos;
//   else
//     return static_cast<Index>(i);
// }


/**
  Constant index operator. We redifine this here so that we can have
  range checking by assert.
    
  \param[in] n Index  
*/
template<class charT>
inline char my_basic_string<charT>::operator[](Index n) const
{
  assert(0<=n);
  assert(n<nelem());
  return basic_string<charT>::operator[](n);
}

/**
  Non-constant index operator. We redifine this here so that we can
  have range checking by assert.
    
  \param[in] n Index
*/
template<class charT>
inline char& my_basic_string<charT>::operator[](Index n)
{
  assert(0<=n);
  assert(n<nelem());
  return basic_string<charT>::operator[](n);
}

// Non-member functions:

// /** Output operator. */
// inline ostream& operator<<(ostream& os, const my_basic_string& v)
// {
//   my_basic_string<base>::const_iterator         i = v.begin();
//   const my_basic_string<base>::const_iterator end = v.end();

//   if ( i!=end )
//     {
//       os << *i;
//       ++i;
//     }

//   for ( ; i!=end; ++i )
//     {
//       os << "\n" << setw(3) << *i;
//     }


//   // Just use the operator of string.
//   operator<<(os,v);
//   return os;
// }


/** The String type for ARTS. Implementation see documentation of
    class my_basic_string. */
typedef my_basic_string<char> String;

// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Strings. */
typedef Array<String> ArrayOfString;

//
// We don't use this function, we use string.find() instead.
//
// //! Find first occurance.
// /*!
//   This returns the index of the first occurance of w in 
//   string x.  
// 
//   A return value of -1 indicates that no matching element was found.
// 
//   \return   The index of the thing we looked for.
//   \param  x The string to search.
//   \param w  The character to look for.
// 
//   \author Stefan Buehler
//   \date   2002-12-06
// */
// template <class base>
// Index find_first( const my_basic_string<base>& x,
//                   const base& w )
// {
//   for ( Index i=0; i<x.nelem(); ++i )
//     if ( w == x[i] )
//       return i;
// 
//   return -1;
// }


#endif  // string_h

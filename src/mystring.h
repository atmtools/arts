/** 
    \file   array.h

   This file contains the definition of String, the ARTS string class.

   \author Stefan Buehler
   \date   2001-09-14
*/

#ifndef string_h
#define string_h

#include <string>
#include "arts.h"

/**
   The implementation for String, the ARTS string class. 

   This adds some additional functionality to the standard stl string
   class, notably:

   a) Range checking by assert

   b) nelem() member function, return the size of the String of type
   Index. 

   The type std::string is just a typedef for
   std::basic_string<char>. Therefore, to make everything work
   correctly, we have to derive our own class from std::basic_string,
   not from std::string directly.
*/
template<class charT>
class my_basic_string : public std::basic_string<charT>
{
public:
  // Constructors:
  my_basic_string();
  explicit my_basic_string(Index n, char c=' ');
  my_basic_string(const basic_string<charT>& A,
		  Index pos=0,
		  Index npos=my_basic_string<charT>::npos);
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
  const char operator[](Index n) const;
  char& operator[](Index n);

  /** Define npos: */
  static const Index npos = static_cast<Index>(std::basic_string<charT>::npos);
};


// Member functions for my_basic_string:


// Constructors:

/** Default constructor. */
template<class charT>
inline my_basic_string<charT>::my_basic_string() : std::basic_string<charT>()  
{ /* Nothing to do here. */ };

/** Constructor setting size. You may give as a second argument a
    character with which to fill the new string. Per default this is
    zero. 
    
    \param n Number of characters
    \param c Optional fill character
*/
template<class charT>
inline my_basic_string<charT>::my_basic_string(Index n, char c) :
  std::basic_string<charT>(n,c) 
{ /* Nothing to do here. */ };

/** Construnctor from another my_basic_string. */
// template<class charT>
// inline my_basic_string<charT>::my_basic_string(const my_basic_string& A) : std::basic_string<charT>(A) 
// { /* Nothing to do here. */ };

/** Construnctor from a basic_string. This is important for handling
    of expressions like this to work correctly:

    String a = b+'.'+c 

    As for std::basic_string, this constructor can also be used to
    initialize the new string from a subrange of the original string. 

    \param A The original string
    \param pos Start position (0 means from the beginning)
    \param npos How many characters to copy
*/
template<class charT>
inline my_basic_string<charT>::my_basic_string(const basic_string<charT>& A,
					       Index pos,
					       Index npos)
{ 
  // Range checks:
  assert(0<=pos);		// Start index must be 0 or greater 0.

//   cout << "A = " << A << "\n";
//   cout << "pos = " << pos << "\n";
//   cout << "size = " << A.size() << "\n";

  assert(static_cast<std::basic_string<charT>::size_type>(pos)<A.size());	
  // At most the last element of the original string.

  assert( npos==my_basic_string<charT>::npos ||
	  ( (npos >= 0) &&
	    (static_cast<std::basic_string<charT>::size_type>(npos)<=(A.size()-pos))
	    )
	  );  // Number of characters to copy must be at the most the
	      // number left. -1 means all remaining characters. 

  // The assertions look complicated, because we have to cast pos and
  // npos to the unsigned size type of basic string to avoid warning
  // messages from the compiler. Both casts are save, because previous
  // assertions check that pos and npos are positive. (The allowed
  // case npos -1 (=my_basic_string<charT>::npos) is also handled
  // correctly.)

  std::basic_string<charT>::operator=(std::basic_string<charT>(A,pos,npos));

};

/** Constructor from a C-style char array. */
template<class charT>
inline my_basic_string<charT>::my_basic_string(const char A[]) : std::basic_string<charT>(A) 
{ /* Nothing to do here. */ };


/** Assignment from another my_basic_string.

    Note that my_basic_string behaves differently from Vector, Matrix, and
    Array here. The two partners do not have to have the same
    size. Size of the target string is adjusted automatically, just as
    std::string does it.
*/
template<class charT>
inline my_basic_string<charT>& my_basic_string<charT>::operator=(const my_basic_string<charT>& A)
{
  std::basic_string<charT>::operator=(A);
  return *this;
}

/** Number of elements. */
template<class charT>
inline Index my_basic_string<charT>::nelem() const
{ 
  size_t s = size();
  assert(s<LONG_MAX);
  return static_cast<long>(s);
}

// /** Find function for char.

//     \param c What character to find.
//     \return Position of c, or npos if not found.

//     Unfortunately, the std::basid_string.find() functions returns npos
//     when the character is not found. This is -1, but assigned to a
//     positive type! Gives a very high number. If you add 1 to this
//     number you get zero again, so it does indeed bahave like -1 in a
//     way. Yuck! With Index this does not work. */
// template<class charT>
// inline Index my_basic_string<charT>::find(char c)
// {
//   std::basic_string<charT>::size_type i = std::basic_string<charT>::find(c);
//   if ( i == std::basic_string<charT>::npos )
//     return npos;
//   else
//     return static_cast<Index>(i);
// }

// /** Find function for string.

//     \param c What string to find.
//     \return Position of c, or npos if not found.

//     Unfortunately, the std::basid_string.find() functions returns npos
//     when the character is not found. This is -1, but assigned to a
//     positive type! Gives a very high number. If you add 1 to this
//     number you get zero again, so it does indeed bahave like -1 in a
//     way. Yuck! With Index this does not work. */
// template<class charT>
// inline Index my_basic_string<charT>::find(const my_basic_string<charT>& c)
// {
//   std::basic_string<charT>::size_type i = std::basic_string<charT>::find(c);
//   if ( i == std::basic_string<charT>::npos )
//     return npos;
//   else
//     return static_cast<Index>(i);
// }


/** Constant index operator. We redifine this here so that we can have
    range checking by assert. */
template<class charT>
inline const char my_basic_string<charT>::operator[](Index n) const
{
  assert(0<=n);
  assert(n<nelem());
  return std::basic_string<charT>::operator[](n);
}

/** Non-constant index operator. We redifine this here so that we can
    have range checking by assert. */
template<class charT>
inline char& my_basic_string<charT>::operator[](Index n)
{
  assert(0<=n);
  assert(n<nelem());
  return std::basic_string<charT>::operator[](n);
}

// Non-member functions:

// /** Output operator. */
// inline std::ostream& operator<<(std::ostream& os, const my_basic_string& v)
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


//   // Just use the operator of std::string.
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


#endif  // string_h

/* Copyright (C) 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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
    \file   array.h

   This file contains the definition of Array.

   \author Stefan Buehler
   \date   2001-09-05
*/

#ifndef array_h
#define array_h

#include <vector>
#include <iomanip>
#include <iostream>
#include <climits>
#include "arts.h"

/** This can be used to make arrays out of anything. 

    A simple #define does not do for this, since I have to implement
    member functions like nelem, in order to be consistent with
    Vector.

    Because constructors are not inherited, I have to re-define all
    costructors. In addition to the constructors here, explicit
    constructors are provided by the derived class MakeArray.
*/
template<class base>
class Array : public std::vector<base>
{
public:
  // Constructors:
  Array()                     : std::vector<base>()  { /* Nothing to do here. */ };
  explicit Array(Index n)     : std::vector<base>(n) { /* Nothing to do here. */ };
  Array(Index n, const base& fill);
  Array(const Array<base>& A) : std::vector<base>(A) { /* Nothing to do here. */ };

  // Assignment operators:
  Array& operator=(base x);
  Array& operator=(const Array<base>& A);

  // Number of elements:
  Index nelem() const;

  // Index operators:
  const base& operator[](Index n) const;
  base& operator[](Index n);
};



// Member functions for Array:

/** Constructor filling with constant value. */
template<class base>
inline Array<base>::Array(Index n, const base& fill) :
  std::vector<base>(n)
{
  // Use std::fill to fill.
  std::fill(this->begin(),this->end(),fill);
};


/** Assignment from base type (fill entire Array with this value). */
template<class base>
inline Array<base>& Array<base>::operator=(base x) 
{
  std::fill(this->begin(),this->end(),x);
  return *this;
}

/** Assignment from another Array. 

    The Behavior of this one is a bit special: If the size of the
    target Array is 0 then it will be automatically resized to match
    (this is needed to have the correct initialization for constructed
    classes that use the assignment operator to initialize their
    data). 

    In all other cases sizes must match exactly!
*/
template<class base>
inline Array<base>& Array<base>::operator=(const Array<base>& A)
{
  //  cout << "size this / A = " << size() << " / " << A.size() << "\n";
  if ( 0==this->size() )
    this->resize(A.size()); // Adjust if previously empty.
  else
    assert( this->size()==A.size() );    // Otherwise check that sizes are compatible.

  std::copy( A.begin(), A.end(), this->begin() );
  return *this;
}

/** Number of elements. */
template<class base>
inline Index Array<base>::nelem() const
{ 
  size_t s = this->size();
  assert(s<LONG_MAX);
  return static_cast<INDEX>(s);
}

/** Constant index operator. We redifine this here so that we can have
    range checking by assert. */
template<class base>
inline const base& Array<base>::operator[](Index n) const
{
  assert(0<=n);
  assert(n<nelem());
  return std::vector<base>::operator[](n);
}

/** Non-constant index operator. We redifine this here so that we can
    have range checking by assert. */
template<class base>
inline base& Array<base>::operator[](Index n)
{
  assert(0<=n);
  assert(n<nelem());
  return std::vector<base>::operator[](n);
}

// Non-member functions:

/** Output operator. */
template<class base>
inline std::ostream& operator<<(std::ostream& os, const Array<base>& v)
{
  typename Array<base>::const_iterator         i = v.begin();
  const typename Array<base>::const_iterator end = v.end();

  if ( i!=end )
    {
      os << setw(3) << *i;
      ++i;
    }

  for ( ; i!=end; ++i )
    {
      os << "\n" << setw(3) << *i;
    }

  return os;
}

/** Max function. */
template<class base>
inline base max(const Array<base>& x)
{ 
  // Initial value for max:
  base max = x[0];

  typename Array<base>::const_iterator       xi = x.begin();
  const typename Array<base>::const_iterator xe = x.end();

  for ( ; xi!=xe ; ++xi )
    {
      if ( *xi > max )
        max = *xi;
    }

  return max;
}

/** Min function. */
template<class base>
inline base min(const Array<base>& x)
{ 
  // Initial value for min:
  base min = x[0];

  typename Array<base>::const_iterator       xi = x.begin();
  const typename Array<base>::const_iterator xe = x.end();

  for ( ; xi!=xe ; ++xi )
    {
      if ( *xi < min )
        min = *xi;
    }

  return min;
}


// It is not a good idea to put all the predefined array types in one
// place. If I do this than a file cannot use one without defining all
// the others.

#endif  // array_h

/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   complex.h
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Saturday, December 14 2002
  
  \brief  A class implementing complex numbers for ARTS.
*/

#ifndef complex_h
#define complex_h

#include <iostream>
#include "arts.h"

//! Complex numbers.
/*! 
  The complex numbers implemented by this class should behave in the
  normal way. 

  IMPORTANT: In analogy to normal C++ float and double, the default
  constructor does not initialize the new variable!
*/
class Complex {
public:

  // Constructors:
  // -------------

  //! Default constructor.
  /*!
    IMPORTANT: In analogy to normal C++ float and double, the
    default constructor does not initialize the new variable!
  */
  Complex(){ /* Nothing to do. */ }

  //! Constructor from explicit real and imaginary part.
  /*! 
  \param re Real part.
  \param im Imaginary part.
  */
  Complex( Numeric re, Numeric im ) : mre(re), mim(im) { /* Nothing to do. */ }

  //! Constructor from real part only.
  /*! 
    This sets the imaginary part to zero.
    
    \param re Real part.
  */
  Complex( Numeric re ) : mre(re), mim(0) { /* Nothing to do. */ }

  //! Copy constructor.
  /*!     
  \param x Another complex number.
  */
  Complex( const Complex& x ) : mre(x.mre), mim(x.mim) { /* Nothing to do. */ }

  //! Real part.
  Numeric Re() const
  {
    return mre;
  }

  //! Real part as lvalue.
  Numeric& Re()
  {
    return mre;
  }

  //! Imaginary part.
  Numeric Im() const
  {
    return mim;
  }

  //! Imaginary part as lvalue.
  Numeric& Im()
  {
    return mim;
  }

  //! Absolute value.
  Numeric Abs() const
  {
    return sqrt( mre*mre + mim*mim );
  }

  //! Phase.
  /*! 
    \return The phase.

    \author FIXME
    \date   FIXME
  */
  Numeric Pha() const
  {
    // FIXME: Nikolay, implement this, please. Don't forget to put a
    // test in test_complex.cc to check whether it works correctly.

    cout << "Not implemented.";
    exit(1);
  }

  //! Addition of another complex number.
  /*! 
    \param x Another complex number to add to this one.
    
    \return We also return *this, to be consistent with other C += operators. 

    \author Stefan Buehler
    \date   2002-12-16
  */
  Complex operator+=(const Complex& x)
  {
    mre += x.mre;
    mim += x.mim;
    return *this;
  }

  /*  FIXME: Nikolay, please implement -=, *=, and /= in similar
      fashion. Don't forget \author and \date tags. Don't forget to
      put tests in test_complex.cc to see if your operators work
      correctly. In the case of /= put in an assertion to detect
      division by zero!
  */

  //! Addition of too complex numbers.
  /*!
    This uses the += operator to add two complex numbers.
  
    \param x The other number to add to *this.
    
    \return The result.
    
    \author Stefan Buehler
    \date   2002-12-16
  */
  Complex operator+( const Complex& x )
  {
    // Create a dummy and initialize with *this:
    Complex a(*this);
    
    // Add x:
    a+=x;

    // return result:
    return a;
  }
  
  /*  FIXME: Nikolay, please implement -, *, and / in similar
      fashion. Note that you should just use the -=, *=, and /=
      operators to achieve this, as I have done in my example.

      Don't forget \author and \date tags. Don't forget to
      put tests in test_complex.cc to see if your operators work
      correctly. 

      When you are done, please check the Doxygen documentation
      generated for file complex.cc, after you said "make" in the
      doc/doxygen directory.
  */

private:
  //! Real part.
  Numeric mre;
  //! Imaginary part.
  Numeric mim;
};


std::ostream& operator<<( std::ostream& os,
			  const Complex& x );

#endif // complex_h

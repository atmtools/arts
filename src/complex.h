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
#include "cmath"
#include "assert.h"
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

  //! Real part as value.
  Numeric& Re()
  {
    return mre;
  }

  //! Imaginary part.
  Numeric Im() const
  {
    return mim;
  }

  //! Imaginary part as value.
  Numeric& Im()
  {
    return mim;
  }
  
  //asolute value
  //  Numeric Abs()
  //{
  // return sqrt( (mre*mre)+(mim*mim))
  //  }
  //asolute value as value
  // Numeric& Abs()
  //{
   //return sqrt( (mre*mre)+(mim*mim))
   //}

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

  //! Subtraction of another complex number.
  /*! 
    \param x Another complex number to add to this one.
    
    \return We also return *this, to be consistent with other C += operators. 

    \author Nikolay Koulev
    \date   2002-12-16
  */
  Complex operator-=(const Complex& x)
  {
    mre -= x.mre;
    mim -= x.mim;
    return *this;
  }
  //! Multiplication of another complex number.
  /*! 
    \param x Another complex number to add to this one.
    
    \return We also return *this, to be consistent with other C += operators. 

    \author Nikolay Koulev, Thomas Kuhn
    \date   2002-12-16
  */
   Complex operator*=(const Complex& x)
  { 
    mre = (mre*x.mre) - (mim*x.mim);
    mim = (mre*x.mim) + (mim*x.mre);
    return *this;
  }
  //! Division of another complex number.
  /*! 
    \param x Another complex number to add to this one.
    
    \return We also return *this, to be consistent with other C += operators. 

    \author Nikolay Koulev, Thomas Kuhn
    \date   2002-12-16
  */
   Complex operator/=(const Complex& x)
  {
    mre = (mre*x.mre+mim*x.mim) / (x.mre*x.mre+x.mim*x.mim);
    mim = (mim*x.mre-mre*x.mim) / (x.mre*x.mre+x.mim*x.mim);
    // assert ((x.mre*x.mre+x.mim*x.mim)!=0)
    return *this;
  }


  /*  FIXME: Nikolay, please implement -=, *=, and /= in similar
      fashion. Don't forget \author and \date tags. Don't forget to
      put tests in test_complex.cc to see if your operators work
      correctly. In the case of /= put in an assertion to detect
      division by zero!
  */

  //! Addition of two complex numbers.
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
    a += x;

    // return result:
    return a;
  }
   //! Subtraction of two complex numbers.
  /*!
    This uses the += operator to add two complex numbers.
  
    \param x The other number to add to *this.
    
    \return The result.
    
    \author Nikolay Koulev
    \date   2002-12-16
  */
  Complex operator-( const Complex& x )
  {
    // Create a dummy and initialize with *this:
    Complex a(*this);
    
    // Subtract x:
    a -= x;

    // return result:
    return a;
  } 
  //! Multiplication of two complex numbers.
  /*!
    This uses the += operator to add two complex numbers.
  
    \param x The other number to add to *this.
    
    \return The result.
    
    \author Nikolay Koulev
    \date   2002-12-16
  */
  Complex operator*( const Complex& x )
  {
    // Create a dummy and initialize with *this:
    Complex a(*this);
    
    // Multiply with x:
    a*=x;

    // return result:
    return a;
   } 
  //! Multiplication of two complex numbers.
  /*!
    This uses the += operator to add two complex numbers.
  
    \param x The other number to add to *this.
    
    \return The result.
    
    \author Nikolay Koulev
    \date   2002-12-16
  */
  Complex operator/( const Complex& x )
  {
    // Create a dummy and initialize with *this:
    Complex a(*this);
    
    // Devide by x:
    a/=x;

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

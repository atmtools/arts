/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   complex.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-12-19
  
  \brief  A class implementing complex numbers for ARTS.
*/

#include "blas.h"
#include "complex.h"
#include "exceptions.h"
#include <cstring>
#include <cmath>

using std::setw;
using std::runtime_error;

std::complex<float> operator+ (const double& d, const std::complex<float>& c)
{
    return (float(d) + c);
}

std::complex<float> operator* (const double& d, const std::complex<float>& c)
{
    return (float(d) * c);
}


std::complex<float> operator+ (const std::complex<float>& c, const double& d)
{
    return (c + float(d));
}

std::complex<float> operator* (const std::complex<float>& c, const double& d)
{
    return (c * float(d));
}



std::complex<double> operator+ (const float& f, const std::complex<double>& c)
{
    return (double(f) + c);
}

std::complex<double> operator* (const float& f, const std::complex<double>& c)
{
    return (double(f) * c);
}

std::complex<double> operator+ (const std::complex<double>& c, const float& f)
{
    return (c + double(f));
}

std::complex<double> operator* (const std::complex<double>& c, const float& f)
{
    return (c * double(f));
}


// Define the global joker object:
extern const Joker joker;

static const Numeric RAD2DEG = 57.295779513082323;

// Functions for ComplexRange:
// --------------------

/* Explicit constructor. 
 * 
 * \param Start must be >= 0.
 * 
 * \param Extent also. Although internally negative extent means "to the end",
 * this can not be created this way, only with the joker. Zero
 * extent is allowed, though, which corresponds to an empty range.
 * 
 * \param Stride can be anything. It can be omitted, in which case the
 * default value is 1. */
ComplexRange::ComplexRange(Index start, Index extent, Index stride) :
mstart(start), mextent(extent), mstride(stride)
{
    // Start must be >= 0:
    assert( 0<=mstart );
    // Extent also. Although internally negative extent means "to the end",
    // this can not be created this way, only with the joker. Zero
    // extent is allowed, though, which corresponds to an empty range.
    assert( 0<=mextent );
    // Stride can be anything except 0.
    // SAB 2001-09-21: Allow 0 stride.
    //  assert( 0!=mstride);
}

/** Constructor with joker extent. Depending on the sign of stride,
 *   this means "to the end", or "to the beginning". */
ComplexRange::ComplexRange(Index start, Joker, Index stride) :
mstart(start), mextent(-1), mstride(stride)
{
    // Start must be >= 0:
    assert( 0<=mstart );
}

/** Constructor with just a joker. This means, take everything. You
 *   can still optionally give a stride, though. This constructor is
 *   just shorter notation for Range(0,joker) */
ComplexRange::ComplexRange(Joker, Index stride) :
mstart(0), mextent(-1), mstride(stride)
{
    // Nothing to do here.
}

/** Constructor which converts a range with joker to an explicit
 *   range.
 * 
 *   \param max_size The maximum allowed size of the vector. 
 *   \param r The new range, with joker. */
ComplexRange::ComplexRange(Index max_size, const ComplexRange& r) :
mstart(r.mstart),
mextent(r.mextent),
mstride(r.mstride)
{
    // Start must be >= 0:
    assert( 0<=mstart );
    // ... and < max_size:
    assert( mstart<max_size );
    
    // Stride must be != 0:
    assert( 0!=mstride);
    
    // Convert negative extent (joker) to explicit extent
    if ( mextent<0 )
    {
        if ( 0<mstride )
            mextent = 1 + (max_size-1-mstart)/mstride;
        else
            mextent = 1 + (0-mstart)/mstride;
    }
    else
    {
        #ifndef NDEBUG
        // Check that extent is ok:
        Index fin = mstart+(mextent-1)*mstride;
        assert( 0   <= fin );
        assert( fin <  max_size );
        #endif
    }
}


/** Constructor of a new range relative to an old range. The new range
 *   may contain -1 for the stride, which acts as a joker.
 * 
 *   \param p Previous range.
 *   \param n New range. */
ComplexRange::ComplexRange(const ComplexRange& p, const ComplexRange& n) :
mstart(p.mstart + n.mstart*p.mstride),
mextent(n.mextent),
mstride(p.mstride*n.mstride)
{
    // We have to juggle here a bit with previous, new, and resulting
    // quantities. I.e.;
    // p.mstride: Previous stride
    // n.mstride: New stride (as specified)
    // mstride:   Resulting stride (old*new)
    
    // Get the previous final element:
    Index prev_fin = p.mstart+(p.mextent-1)*p.mstride;
    
    // Resulting start must be >= previous start:
    assert( p.mstart<=mstart );
    // and <= prev_fin:
    assert( mstart<=prev_fin );
    
    // Resulting stride must be != 0:
    assert( 0!=mstride);
    
    // Convert negative extent (joker) to explicit extent
    if ( mextent<0 )
    {
        if ( 0<mstride )
            mextent = 1 + (prev_fin-mstart)/mstride;
        else
            mextent = 1 + (p.mstart-mstart)/mstride;
    }
    else
    {
        #ifndef NDEBUG
        // Check that extent is ok:
        Index fin = mstart+(mextent-1)*mstride;
        assert( p.mstart <= fin      );
        assert( fin      <= prev_fin );
        #endif
    }
    
}

// Functions for ConstComplexVectorView:
// ------------------------------

//! Returns true if variable size is zero.
bool ConstComplexVectorView::empty() const
{
    return (nelem() == 0);
}

/** Returns the number of elements.  The names `size' and `length'
    are already used by STL functions returning size_t. To avoid
    confusion we choose the name `nelem'. This is also more
    consistent with `nrow' and `ncol' for matrices.
    
    The value range of long, which is used to store the index is on a
    PC from -2147483648 to 2147483647. This means that a 15GB large
    array of float can be addressed with this index. So the extra bit
    that size_t has compared to long is not needed. */
Index ConstComplexVectorView::nelem() const
{
  return mrange.mextent;
}

/** The sum of all elements of a Vector. */
Complex ConstComplexVectorView::sum() const
{
  Complex s=0;
  ConstComplexIterator1D i = begin();
  const ConstComplexIterator1D e = end();

  for ( ; i!=e; ++i )
    s += *i;

  return s;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Vector. This allows
    correct recursive behavior.  */
ConstComplexVectorView ConstComplexVectorView::operator[](const ComplexRange& r) const
{
    return ConstComplexVectorView(mdata, mrange, r);
}

/** Return const iterator to first element. */
ConstComplexIterator1D ConstComplexVectorView::begin() const
{
    return ConstComplexIterator1D(mdata+mrange.mstart, mrange.mstride);
}

/** Return const iterator behind last element. */
ConstComplexIterator1D ConstComplexVectorView::end() const
{
    return ConstComplexIterator1D( mdata +
                          mrange.mstart +
                          (mrange.mextent)*mrange.mstride,
                          mrange.mstride );
}

/** Conversion to const 1 column matrix. */
ConstComplexVectorView::operator ConstComplexMatrixView() const
{
    return ConstComplexMatrixView(mdata,mrange,ComplexRange(0,1));
}

/** A special constructor, which allows to make a ConstComplexVectorView from
    a scalar.

    This one is a bit tricky: We have to cast away the arguments const
    qualifier, because mdata is not const. This should be safe, since
    there are no non-const methods for ConstComplexVectorView.
*/
ConstComplexVectorView::ConstComplexVectorView(const Complex& a) :
  mrange(0,1), mdata(&const_cast<Complex&>(a))
{
  // Nothing to do here.
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Vector. */
ConstComplexVectorView::ConstComplexVectorView() :
  mrange(0,0), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
ConstComplexVectorView::ConstComplexVectorView( Complex *data,
                                         const ComplexRange& range) :
  mrange(range),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub ranges from
    sub ranges. That means that the new range has to be interpreted
    relative to the original range. The new range may contain -1 for
    the extent which acts as a joker. However, the used Range
    constructor converts this to an explicit range, consistent with
    the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
ConstComplexVectorView::ConstComplexVectorView( Complex *data,
                                                const ComplexRange& p,
                                                const ComplexRange& n) :
  mrange(p,n),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the vector. The iterators know which part of the vector
    is `active', and also the stride. */
std::ostream& operator<<(std::ostream& os, const ConstComplexVectorView& v)
{
    ConstComplexIterator1D i=v.begin();
    const ConstComplexIterator1D end=v.end();

  if ( i!=end )
    {
      os << setw(3) << *i;
      ++i;
    }
  for ( ; i!=end ; ++i )
    {
      os << " " << setw(3) << *i;
    }

  return os;
}



// Functions for ComplexVectorView:
// ------------------------

/** Bail out immediately if somebody tries to create a ComplexVectorView from
 a const Complex*Vector. */
ComplexVectorView::ComplexVectorView (const ComplexVector&)
{
    throw runtime_error("Creating a ComplexVectorView from a const ComplexVector is not allowed.");
  // This is not really a runtime error, but I don't want to start
  // producing direct output from inside matpack. And just exiting is
  // not so nice. 
  // If you see this error, there is a bug in the code, not in the
  // ARTS input.
}

/** Create ComplexVectorView from a ComplexVector. */
ComplexVectorView::ComplexVectorView (ComplexVector& v)
{
  mdata = v.mdata;
  mrange = v.mrange;
}

/** Const index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  Has to be redifined here, because the
    one from ConstComplexVectorView is hidden. */
ConstComplexVectorView ComplexVectorView::operator[](const ComplexRange& r) const
{
    return ConstComplexVectorView::operator[](r);
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  */
ComplexVectorView ComplexVectorView::operator[](const ComplexRange& r)
{
    return ComplexVectorView(mdata, mrange, r);
}

/** Return const iterator to first element. Has to be redefined here,
    since it is hiden by the non-const operator of the derived
    class.*/
ConstComplexIterator1D ComplexVectorView::begin() const
{
    return ConstComplexVectorView::begin();
}

/** Return const iterator behind last element. Has to be redefined
    here, since it is hiden by the non-const operator of the derived
    class.*/
ConstComplexIterator1D ComplexVectorView::end() const
{
    return ConstComplexVectorView::end();
}

/** Return iterator to first element. */
ComplexIterator1D ComplexVectorView::begin()
{
    return ComplexIterator1D(mdata+mrange.mstart, mrange.mstride);
}

/** Return iterator behind last element. */
ComplexIterator1D ComplexVectorView::end()
{
    return ComplexIterator1D( mdata +
                     mrange.mstart +
                     (mrange.mextent)*mrange.mstride,
                     mrange.mstride );
}

/** Assignment operator. This copies the data from another VectorView
 t o this Complex*VectorView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this VectorView by
    setting its range. */
ComplexVectorView& ComplexVectorView::operator=(const ConstComplexVectorView& v)
{
  //  cout << "Assigning VectorView from ConstVectorView.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);

  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assignment from ComplexVectorView to ComplexVectorView. This is a tricky
 one. The problem is that since ComplexVectorView is derived from
 ConstComplexVectorView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
ComplexVectorView& ComplexVectorView::operator=(const ComplexVectorView& v)
{
    //  cout << "Assigning ComplexVectorView from ComplexVectorView.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);

  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assignment from ComplexVector. This is important to avoid a bug when
 assigning a Vector to a Complex*VectorView. */
ComplexVectorView& ComplexVectorView::operator=(const ComplexVector& v)
{
    //  cout << "Assigning ComplexVectorView from Vector.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);
  
  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assigning a scalar to a ComplexVectorView will set all elements to this
    value. */
ComplexVectorView& ComplexVectorView::operator=(Complex x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
ComplexVectorView ComplexVectorView::operator*=(Complex x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
    *i *= x;
  return *this;
}

/** Multiplication by scalar. */
ComplexVectorView ComplexVectorView::operator*=(Numeric x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
        *i *= x;
    return *this;
}

/** Division by scalar. */
ComplexVectorView ComplexVectorView::operator/=(Complex x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
    *i /= x;
  return *this;
}

/** Division by scalar. */
ComplexVectorView ComplexVectorView::operator/=(Numeric x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
        *i /= x;
    return *this;
}

/** Addition of scalar. */
ComplexVectorView ComplexVectorView::operator+=(Complex x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
    *i += x;
  return *this;
}

/** Addition of scalar. */
ComplexVectorView ComplexVectorView::operator+=(Numeric x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
        *i += x;
    return *this;
}

/** Subtraction of scalar. */
ComplexVectorView ComplexVectorView::operator-=(Complex x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
    *i -= x;
  return *this;
}

/** Subtraction of scalar. */
ComplexVectorView ComplexVectorView::operator-=(Numeric x)
{
    const ComplexIterator1D e=end();
    for ( ComplexIterator1D i=begin(); i!=e ; ++i )
        *i -= x;
    return *this;
}

/** Element-vise multiplication by another vector. */
ComplexVectorView ComplexVectorView::operator*=(const ConstComplexVectorView& x)
{
    assert( nelem()==x.nelem() );
    
    ConstComplexIterator1D s=x.begin();
    
    ComplexIterator1D i=begin();
    const ComplexIterator1D e=end();
    
    for ( ; i!=e ; ++i,++s )
        *i *= *s;
    return *this;
}

/** Element-vise multiplication by another vector. */
ComplexVectorView ComplexVectorView::operator*=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  ComplexIterator1D i=begin();
  const ComplexIterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i *= *s;
  return *this;
}

/** Element-vise division by another vector. */
ComplexVectorView ComplexVectorView::operator/=(const ConstComplexVectorView& x)
{
    assert( nelem()==x.nelem() );
    
    ConstComplexIterator1D s=x.begin();
    
    ComplexIterator1D i=begin();
    const ComplexIterator1D e=end();
    
    for ( ; i!=e ; ++i,++s )
        *i /= *s;
    return *this;
}

/** Element-vise division by another vector. */
ComplexVectorView ComplexVectorView::operator/=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  ComplexIterator1D i=begin();
  const ComplexIterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i /= *s;
  return *this;
}

/** Element-vise addition of another vector. */
ComplexVectorView ComplexVectorView::operator+=(const ConstComplexVectorView& x)
{
    assert( nelem()==x.nelem() );
    
    ConstComplexIterator1D s=x.begin();
    
    ComplexIterator1D i=begin();
    const ComplexIterator1D e=end();
    
    for ( ; i!=e ; ++i,++s )
        *i += *s;
    return *this;
}

/** Element-vise addition of another vector. */
ComplexVectorView ComplexVectorView::operator+=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  ComplexIterator1D i=begin();
  const ComplexIterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i += *s;
  return *this;
}

/** Element-vise subtraction of another vector. */
ComplexVectorView ComplexVectorView::operator-=(const ConstComplexVectorView& x)
{
    assert( nelem()==x.nelem() );
    
    ConstComplexIterator1D s=x.begin();
    
    ComplexIterator1D i=begin();
    const ComplexIterator1D e=end();
    
    for ( ; i!=e ; ++i,++s )
        *i -= *s;
    return *this;
}

/** Element-vise subtraction of another vector. */
ComplexVectorView ComplexVectorView::operator-=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  ComplexIterator1D i=begin();
  const ComplexIterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i -= *s;
  return *this;
}

/** Conversion to 1 column matrix. */
ComplexVectorView::operator ComplexMatrixView()
{
    // The old version (before 2013-01-18) of this was:
    //    return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
    // Bus this was a bug! The problem is that the matrix index operator adds
    // the mstart from both row and columm range object to mdata

    return ComplexMatrixView(mdata,mrange,ComplexRange(0,1));
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  VectorView is not pointing to the beginning of a Vector or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Complex * ComplexVectorView::get_c_array() const
{
    if (mrange.mstart != 0 || mrange.mstride != 1)
        throw runtime_error("A ComplexVectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");

    return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  VectorView is not pointing to the beginning of a Vector or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Complex * ComplexVectorView::get_c_array()
{
    if (mrange.mstart != 0 || mrange.mstride != 1)
        throw runtime_error("A VectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");
    
  return mdata;
}

/** A special constructor, which allows to make a VectorView from
    a scalar.
*/
ComplexVectorView::ComplexVectorView(Complex& a) :
ConstComplexVectorView(a)
{
  // Nothing to do here.
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Vector. */
ComplexVectorView::ComplexVectorView() :
ConstComplexVectorView()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
ComplexVectorView::ComplexVectorView(Complex *data,
                       const ComplexRange& range) :
                       ConstComplexVectorView(data,range)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub ranges from
    sub ranges. That means that the new range has to be interpreted
    relative to the original range. The new range may contain -1 for
    the extent which acts as a joker. However, the used Range
    constructor converts this to an explicit range, consistent with
    the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
ComplexVectorView::ComplexVectorView(Complex *data,
                                     const ComplexRange& p,
                                     const ComplexRange& n) :
                                     ConstComplexVectorView(data,p,n)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors. */
void copy(ConstComplexIterator1D origin,
          const ConstComplexIterator1D& end,
          ComplexIterator1D target)
{
  if (origin.mstride == 1 && target.mstride == 1)
    memcpy((void *)target.mx,
           (void *)origin.mx,
           sizeof(Complex) * (end.mx - origin.mx));
  else
    for ( ; origin!=end ; ++origin,++target )
      *target = *origin;
}

/** Copy a scalar to all elements. */
void copy(Complex x,
          ComplexIterator1D target,
          const ComplexIterator1D& end)
{
  for ( ; target!=end ; ++target )
    *target = x;
}

// Functions for Vector:
// ---------------------

/** Default constructor. */
ComplexVector::ComplexVector() 
{
    // Nothing to do here
}

/** Constructor setting size. */
ComplexVector::ComplexVector(Index n) :
ComplexVectorView( new Complex[n],
                   ComplexRange(0,n))
{
    // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
ComplexVector::ComplexVector(Index n, Complex fill) :
ComplexVectorView( new Complex[n],
                   ComplexRange(0,n))
{
    // Here we can access the raw memory directly, for slightly
    // increased efficiency:
    const Complex *stop = mdata+n;
    for ( Complex *x=mdata; x<stop; ++x )
        *x = fill;
}

/** Constructor setting size and filling with constant value. */
ComplexVector::ComplexVector(Index n, Numeric fill) :
ComplexVectorView( new Complex[n],
                   ComplexRange(0,n))
{
    // Here we can access the raw memory directly, for slightly
    // increased efficiency:
    const Complex *stop = mdata+n;
    for ( Complex *x=mdata; x<stop; ++x )
        *x = fill;
}

/** Constructor filling with values. 
 * 
 *   Examples:
 * 
 *   Vector v(1,5,1);  // 1, 2, 3, 4, 5
 *   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
 *   Vector v(5,5,-1); // 5, 4, 3, 2, 1
 */
ComplexVector::ComplexVector(Complex start, Index extent, Complex stride) :
ComplexVectorView( new Complex[extent],
                   ComplexRange(0,extent))
{
    // Fill with values:
    Complex x = start;
    ComplexIterator1D       i=begin();
    const ComplexIterator1D e=end();
    for ( ; i!=e; ++i )
    {
        *i = x;
        x += stride;
    }
}

/** Constructor filling with values. 
* 
*   Examples:
* 
*   Vector v(1,5,1);  // 1, 2, 3, 4, 5
*   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
*   Vector v(5,5,-1); // 5, 4, 3, 2, 1
*/
ComplexVector::ComplexVector(Numeric start, Index extent, Complex stride) :
ComplexVectorView( new Complex[extent],
                   ComplexRange(0,extent))
{
    // Fill with values:
    Complex x = start;
    ComplexIterator1D       i=begin();
    const ComplexIterator1D e=end();
    for ( ; i!=e; ++i )
    {
        *i = x;
        x += stride;
    }
}

/** Constructor filling with values. 
 * 
 *   Examples:
 * 
 *   Vector v(1,5,1);  // 1, 2, 3, 4, 5
 *   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
 *   Vector v(5,5,-1); // 5, 4, 3, 2, 1
 */
ComplexVector::ComplexVector(Complex start, Index extent, Numeric stride) :
ComplexVectorView( new Complex[extent],
                   ComplexRange(0,extent))
{
    // Fill with values:
    Complex x = start;
    ComplexIterator1D       i=begin();
    const ComplexIterator1D e=end();
    for ( ; i!=e; ++i )
    {
        *i = x;
        x += stride;
    }
}

/** Constructor filling with values. 
 * 
 *   Examples:
 * 
 *   Vector v(1,5,1);  // 1, 2, 3, 4, 5
 *   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
 *   Vector v(5,5,-1); // 5, 4, 3, 2, 1
 */
ComplexVector::ComplexVector(Numeric start, Index extent, Numeric stride) :
ComplexVectorView( new Complex[extent],
                   ComplexRange(0,extent))
{
    // Fill with values:
    Complex x = start;
    ComplexIterator1D       i=begin();
    const ComplexIterator1D e=end();
    for ( ; i!=e; ++i )
    {
        *i = x;
        x += stride;
    }
}

/** Copy constructor from ComplexVectorView. This automatically sets the size
 *   and copies the data. The vector created will have start zero and
 *   stride 1, independent on how these parameters are set for the
 *   original. So, what is copied is the data, not the shape
 *   of the selection. */
ComplexVector::ComplexVector(const ConstComplexVectorView& v) :
ComplexVectorView( new Complex[v.nelem()],
                   ComplexRange(0,v.nelem()))
{
    copy(v.begin(),v.end(),begin());
}

/** Copy constructor from ComplexVector. This is important to override the
 *   automatically generated shallow constructor. We want deep copies!  */
ComplexVector::ComplexVector(const ComplexVector& v) :
ComplexVectorView( new Complex[v.nelem()],
                   ComplexRange(0,v.nelem()))
{
    copy(v.begin(),v.end(),begin());
}

/** Converting constructor from std::vector. */
ComplexVector::ComplexVector(const std::vector<Complex>& v) :
ComplexVectorView( new Complex[v.size()],
                   ComplexRange(0,v.size()))
{
    std::vector<Complex>::const_iterator vec_it_end = v.end();
    ComplexIterator1D this_it = this->begin();
    for (std::vector<Complex>::const_iterator vec_it = v.begin();
         vec_it != vec_it_end;
    ++vec_it, ++this_it)
         *this_it = *vec_it;
}

/** Converting constructor from std::vector. */
ComplexVector::ComplexVector(const std::vector<Numeric>& v) :
ComplexVectorView( new Complex[v.size()],
                   ComplexRange(0,v.size()))
{
    std::vector<Numeric>::const_iterator vec_it_end = v.end();
    ComplexIterator1D this_it = this->begin();
    for (std::vector<Numeric>::const_iterator vec_it = v.begin();
         vec_it != vec_it_end;
    ++vec_it, ++this_it)
         *this_it = *vec_it;
}

//! Assignment from another Vector.
/*! 
 * While dimensions of VectorViews can not be adjusted, dimensions of
 * Vectors *can* be adjusted. Hence, the behavior of the assignment
 * operator is different.
 * 
 * In this case the size of the target is automatically adjusted. This
 * is important, so that structures containing Vectors are copied
 * correctly. 
 * 
 * This is a deviation from the old ARTS paradigm that sizes must match
 * exactly before copying!
 * 
 * Note: It is sufficient to have only this one version of the
 * assignment (Vector = Vector). It implicitly covers the cases
 * Vector=VectorView, etc, because there is a default constructor for
 * Vector from VectorView. (See C++ Primer Plus, page 571ff.)
 * 
 * \param v The other vector to copy to this one.
 * 
 * \return This vector, by tradition.
 * 
 * \author Stefan Buehler
 * \date   2002-12-19
 */
ComplexVector& ComplexVector::operator=(ComplexVector v)
{
    swap(*this, v);
    return *this;
}

//! Assignment operator from Array<Numeric>.
/*!
 * This copies the data from a Array<Numeric> to this VectorView. The
 * size is adjusted automatically.
 * 
 * Array<Numeric> can be useful to collect things in, because there
 * is a .push_back method for it. Then, after collecting we usually
 * have to transfer the content to a Vector. With this assignment
 * operator that's easy.
 * 
 * \param x The array to assign to this.
 * 
 * \return This vector, by tradition.
 * 
 * \author Stefan Buehler
 * \date   2002-12-19
 */
ComplexVector& ComplexVector::operator=(const Array<Complex>& x)
{
    resize( x.nelem() ); 
    ComplexVectorView::operator=(x);
    return *this;
}

/** Assignment operator from scalar. Assignment operators are not
 *   inherited. */  
ComplexVector& ComplexVector::operator=(Complex x)
{
    ComplexVectorView::operator=(x);
    return *this;
}

// /** Assignment operator from VectorView. Assignment operators are not
//     inherited. */  
// inline Vector& Vector::operator=(const VectorView v)
// {
//   cout << "Assigning Vector from Vector View.\n";
//  // Check that sizes are compatible:
//   assert(mrange.mextent==v.mrange.mextent);
//   VectorView::operator=(v);
//   return *this;
// }

/** Resize function. If the size is already correct this function does
 *   nothing. All data is lost after resizing! The new Vector is not
 *   initialized, so it will contain random values.  */
void ComplexVector::resize(Index n)
{
    assert( 0<=n );
    if ( mrange.mextent != n )
    {
        delete[] mdata;
        mdata = new Complex[n];
        mrange.mstart = 0;
        mrange.mextent = n;
        mrange.mstride = 1;
    }
}


/** Swaps two objects. */
void swap(ComplexVector& v1, ComplexVector& v2)
{
    std::swap(v1.mrange, v2.mrange);
    std::swap(v1.mdata, v2.mdata);
}


/** Destructor for ComplexVector. This is important, since Vector uses new to
 *   allocate storage. */
ComplexVector::~ComplexVector()
{
    delete[] mdata;
}

// Functions for ConstMatrixView:
// ------------------------------

//! Returns true if variable size is zero.
bool ConstComplexMatrixView::empty() const
{
    return (nrows() == 0 || ncols() == 0);
}

/** Returns the number of rows. */
Index ConstComplexMatrixView::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
Index ConstComplexMatrixView::ncols() const
{
  return mcr.mextent;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior.  */
ConstComplexMatrixView ConstComplexMatrixView::operator()(const ComplexRange& r,
                                                          const ComplexRange& c) const
{
    return ConstComplexMatrixView(mdata, mrr, mcr, r, c);
}

/** Const index operator returning a column as an object of type
 C onstComplex*VectorView.

    \param r A range of rows.
    \param c Index of selected column */
ConstComplexVectorView ConstComplexMatrixView::operator()(const ComplexRange& r, Index c) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ConstComplexVectorView(mdata + mcr.mstart + c*mcr.mstride,
                         mrr, r);
}

/** Const index operator returning a row as an object of type
 C onstComplex*VectorView.

    \param r Index of selected row.
    \param c Range of columns */
ConstComplexVectorView ConstComplexMatrixView::operator()(Index r, const ComplexRange& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return ConstComplexVectorView(mdata + mrr.mstart + r*mrr.mstride,
                         mcr, c);
}

/** Return const iterator to first row. */
ConstComplexIterator2D ConstComplexMatrixView::begin() const
{
    return ConstComplexIterator2D(ConstComplexVectorView(mdata+mrr.mstart,
                                         mcr),
                         mrr.mstride);
}

/** Return const iterator behind last row. */
ConstComplexIterator2D ConstComplexMatrixView::end() const
{
    return ConstComplexIterator2D( ConstComplexVectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
                                          mcr),
                          mrr.mstride );
}

//! ComplexMatrix diagonal as vector.
/*!
 R eturns a ConstComplex*MatrixView on the diagonal entries of the matrix. For a given
  (n,m) matrix M the diagonal vector v is the vector of length min{n,m} with entries

       v[i] = M(i,i)

  \return The diagonal vector v.
*/
ConstComplexVectorView ConstComplexMatrixView::diagonal() const
{
    Index n = std::min( mrr.mextent, mcr.mextent );
    return ConstComplexVectorView( mdata + mrr.mstart + mcr.mstart,
                                   ComplexRange( 0, n, mrr.mstride + mcr.mstride ) );
}


/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
ConstComplexMatrixView::ConstComplexMatrixView() :
  mrr(0,0,1), mcr(0,0,1), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
ConstComplexMatrixView::ConstComplexMatrixView( Complex *data,
                                                const ComplexRange& rr,
                                                const ComplexRange& cr) :
  mrr(rr),
  mcr(cr),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pr Previous range.
    \param pc Previous range.
    \param nr New Range.
    \param nc New Range.
  */
ConstComplexMatrixView::ConstComplexMatrixView( Complex *data,
                                                const ComplexRange& pr, const ComplexRange& pc,
                                                const ComplexRange& nr, const ComplexRange& nc) :
  mrr(pr,nr),
  mcr(pc,nc),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the matrix. The iterators know which part of the matrix
    is `active', and also the strides in both directions. This
    function is a bit more complicated than necessary to illustrate
    the concept, because the formating should look nice. This means
    that the first row, and the first element in each row, have to be
    treated individually. */
std::ostream& operator<<(std::ostream& os, const ConstComplexMatrixView& v)
{
  // Row iterators:
    ConstComplexIterator2D ir=v.begin();
    const ConstComplexIterator2D end_row=v.end();

  if ( ir!=end_row )
    {
        ConstComplexIterator1D ic =  ir->begin();
        const ConstComplexIterator1D end_col = ir->end();

      if ( ic!=end_col )
        {
          os << setw(3) << *ic;
          ++ic;
        }
      for ( ; ic!=end_col ; ++ic )
        {
          os << " " << setw(3) << *ic;
        }
      ++ir;
    }
  for ( ; ir!=end_row ; ++ir )
    {
        ConstComplexIterator1D ic =  ir->begin();
        const ConstComplexIterator1D end_col = ir->end();

      os << "\n";
      if ( ic!=end_col )
        {
          os << setw(3) << *ic;
          ++ic;
        }
      for ( ; ic!=end_col ; ++ic )
        {
          os << " " << setw(3) << *ic;
        }
    }

  return os;
}


// Functions for ComplexMatrixView:
// -------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a ComplexMatrix. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
ConstComplexMatrixView ComplexMatrixView::operator()(const ComplexRange& r, const ComplexRange& c) const
{
    return ConstComplexMatrixView::operator()(r,c);  
}

/** Const index operator returning a column as an object of type
    ConstComplexVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r A range of rows.
    \param c Index of selected column */
ConstComplexVectorView ComplexMatrixView::operator()(const ComplexRange& r, Index c) const
{
    return ConstComplexMatrixView::operator()(r,c);
}

/** Const index operator returning a row as an object of type
    ConstComplexVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r Index of selected row.
    \param c Range of columns */
ConstComplexVectorView ComplexMatrixView::operator()(Index r, const ComplexRange& c) const
{
    return ConstComplexMatrixView::operator()(r,c);
}

/** Index operator for subrange. We have to also account for the case,
 t hat *this is already a subrange of a Complex*Matrix. This allows correct
    recursive behavior.  */
ComplexMatrixView ComplexMatrixView::operator()(const ComplexRange& r, const ComplexRange& c)
{
    return ComplexMatrixView(mdata, mrr, mcr, r, c);
}

/** Index operator returning a column as an object of type ComplexVectorView.

    \param r A range of rows.
    \param c Index of selected column */
ComplexVectorView ComplexMatrixView::operator()(const ComplexRange& r, Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ComplexVectorView(mdata + mcr.mstart + c*mcr.mstride,
                    mrr, r);
}

/** Index operator returning a row as an object of type ComplexVectorView.

    \param r Index of selected row.
    \param c Range of columns */
ComplexVectorView ComplexMatrixView::operator()(Index r, const ComplexRange& c)
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r <  mrr.mextent );

  return ComplexVectorView(mdata + mrr.mstart + r*mrr.mstride,
                    mcr, c);
}

/** Return const iterator to first row. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
ConstComplexIterator2D ComplexMatrixView::begin() const
{
    return ConstComplexMatrixView::begin();
}

/** Return const iterator behind last row. */
ConstComplexIterator2D ComplexMatrixView::end() const
{
    return ConstComplexMatrixView::end();
}

/** Return iterator to first row. */
ComplexIterator2D ComplexMatrixView::begin()
{
    return ComplexIterator2D(ComplexVectorView(mdata+mrr.mstart, mcr),
                    mrr.mstride);
}

/** Return iterator behind last row. */
ComplexIterator2D ComplexMatrixView::end()
{
    return ComplexIterator2D( ComplexVectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
                                mcr),
                     mrr.mstride );
}

/** Assignment operator. This copies the data from another ComplexMatrixView
    to this ComplexMatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this ComplexMatrixView by
    setting its range. */
ComplexMatrixView& ComplexMatrixView::operator=(const ConstComplexMatrixView& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from ComplexMatrixView to ComplexMatrixView. This is a tricky
    one. The problem is that since ComplexMatrixView is derived from
    ConstMatrixView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
ComplexMatrixView& ComplexMatrixView::operator=(const ComplexMatrixView& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a ComplexMatrix. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
ComplexMatrixView& ComplexMatrixView::operator=(const ComplexMatrix& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a vector. This copies the data from a ComplexVectorView
    to this ComplexMatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this ComplexMatrixView by
    setting its range. */
ComplexMatrixView& ComplexMatrixView::operator=(const ConstComplexVectorView& v)
{
  // Check that sizes are compatible:
  assert( mrr.mextent==v.nelem() );
  assert( mcr.mextent==1         );
  //  dummy = ConstComplexMatrixView(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstComplexMatrixView dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Assigning a scalar to a MatrixView will set all elements to this
    value. */
ComplexMatrixView& ComplexMatrixView::operator=(Complex x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
ComplexMatrixView& ComplexMatrixView::operator*=(Complex x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
            *c *= x;
    }
    return *this;
}

/** Multiplication by scalar. */
ComplexMatrixView& ComplexMatrixView::operator*=(Numeric x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
        *c *= x;
    }
  return *this;
}

/** Division by scalar. */
ComplexMatrixView& ComplexMatrixView::operator/=(Complex x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
            *c /= x;
    }
    return *this;
}

/** Division by scalar. */
ComplexMatrixView& ComplexMatrixView::operator/=(Numeric x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
        *c /= x;
    }
  return *this;
}

/** Addition of scalar. */
ComplexMatrixView& ComplexMatrixView::operator+=(Complex x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
            *c += x;
    }
    return *this;
}

/** Addition of scalar. */
ComplexMatrixView& ComplexMatrixView::operator+=(Numeric x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
        *c += x;
    }
  return *this;
}

/** Subtraction of scalar. */
ComplexMatrixView& ComplexMatrixView::operator-=(Complex x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
            *c -= x;
    }
    return *this;
}

/** Subtraction of scalar. */
ComplexMatrixView& ComplexMatrixView::operator-=(Numeric x)
{
    const ComplexIterator2D er=end();
    for ( ComplexIterator2D r=begin(); r!=er ; ++r )
    {
        const ComplexIterator1D ec = r->end();
        for ( ComplexIterator1D c = r->begin(); c!=ec ; ++c )
        *c -= x;
    }
  return *this;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Complex *ComplexMatrixView::get_c_array() const
{
    if (mrr.mstart != 0 || mrr.mstride != mcr.mextent
        || mcr.mstart != 0 || mcr.mstride != 1)
        throw std::runtime_error("A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");

    return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Complex *ComplexMatrixView::get_c_array()
{
    if (mrr.mstart != 0 || mrr.mstride != mcr.mextent
        || mcr.mstart != 0 || mcr.mstride != 1)
        throw std::runtime_error("A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");

    return mdata;
}

/** Element-vise multiplication by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator*=(const ConstComplexMatrixView& x)
{
    assert(nrows()==x.nrows());
    assert(ncols()==x.ncols());
    ConstComplexIterator2D  sr = x.begin();
    ComplexIterator2D        r = begin();
    const ComplexIterator2D er = end();
    for ( ; r!=er ; ++r,++sr )
    {
        ConstComplexIterator1D  sc = sr->begin(); 
        ComplexIterator1D        c = r->begin();
        const ComplexIterator1D ec = r->end();
        for ( ; c!=ec ; ++c,++sc )
            *c *= *sc;
    }
    return *this;
}

/** Element-vise multiplication by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator*=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  ComplexIterator2D        r = begin();
  const ComplexIterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      ComplexIterator1D        c = r->begin();
      const ComplexIterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c *= *sc;
    }
  return *this;
}

/** Element-vise division by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator/=(const ConstComplexMatrixView& x)
{
    assert(nrows()==x.nrows());
    assert(ncols()==x.ncols());
    ConstComplexIterator2D  sr = x.begin();
    ComplexIterator2D        r = begin();
    const ComplexIterator2D er = end();
    for ( ; r!=er ; ++r,++sr )
    {
        ConstComplexIterator1D  sc = sr->begin(); 
        ComplexIterator1D        c = r->begin();
        const ComplexIterator1D ec = r->end();
        for ( ; c!=ec ; ++c,++sc )
            *c /= *sc;
    }
    return *this;
}

/** Element-vise division by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator/=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  ComplexIterator2D        r = begin();
  const ComplexIterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      ComplexIterator1D        c = r->begin();
      const ComplexIterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c /= *sc;
    }
  return *this;
}

/** Element-vise addition of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator+=(const ConstComplexMatrixView& x)
{
    assert(nrows()==x.nrows());
    assert(ncols()==x.ncols());
    ConstComplexIterator2D  sr = x.begin();
    ComplexIterator2D        r = begin();
    const ComplexIterator2D er = end();
    for ( ; r!=er ; ++r,++sr )
    {
        ConstComplexIterator1D  sc = sr->begin(); 
        ComplexIterator1D        c = r->begin();
        const ComplexIterator1D ec = r->end();
        for ( ; c!=ec ; ++c,++sc )
            *c += *sc;
    }
    return *this;
}

/** Element-vise addition of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator+=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  ComplexIterator2D        r = begin();
  const ComplexIterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      ComplexIterator1D        c = r->begin();
      const ComplexIterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c += *sc;
    }
  return *this;
}

/** Element-vise subtraction of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator-=(const ConstComplexMatrixView& x)
{
    assert(nrows()==x.nrows());
    assert(ncols()==x.ncols());
    ConstComplexIterator2D  sr = x.begin();
    ComplexIterator2D        r = begin();
    const ComplexIterator2D er = end();
    for ( ; r!=er ; ++r,++sr )
    {
        ConstComplexIterator1D  sc = sr->begin(); 
        ComplexIterator1D        c = r->begin();
        const ComplexIterator1D ec = r->end();
        for ( ; c!=ec ; ++c,++sc )
            *c -= *sc;
    }
    return *this;
}

/** Element-vise subtraction of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator-=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  ComplexIterator2D        r = begin();
  const ComplexIterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      ComplexIterator1D        c = r->begin();
      const ComplexIterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c -= *sc;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Matrix. */
ComplexMatrixView::ComplexMatrixView() :
ConstComplexMatrixView()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by ComplexMatrix to initialize its
    own ComplexMatrixView part. The row range rr must have a
    stride to account for the length of one row. */
ComplexMatrixView::ComplexMatrixView(Complex *data,
                                     const ComplexRange& rr,
                                     const ComplexRange& cr) :
                                     ConstComplexMatrixView(data, rr, cr)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pr Previous range.
    \param pc Previous range.
    \param nr New Range.
    \param nc New Range.
  */
ComplexMatrixView::ComplexMatrixView(Complex *data,
                                     const ComplexRange& pr, const ComplexRange& pc,
                                     const ComplexRange& nr, const ComplexRange& nc) :
                                     ConstComplexMatrixView(data,pr,pc,nr,nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors.

    Origin, end, and target are 2D iterators, marking rows in a
    matrix. For each row the 1D iterator is obtained and used to copy
    the elements. 
*/
void copy(ConstComplexIterator2D origin,
          const ConstComplexIterator2D& end,
          ComplexIterator2D target)
{
  for ( ; origin!=end ; ++origin,++target )
    {
        ConstComplexIterator1D       o = origin->begin();
        const ConstComplexIterator1D e = origin->end();
        ComplexIterator1D            t = target->begin();
      for ( ; o!=e ; ++o,++t )
        *t = *o;
    }
}

/** Copy a scalar to all elements. */
void copy(Complex x,
          ComplexIterator2D target,
          const ComplexIterator2D& end)
{
    for ( ; target!=end ; ++target )
    {
        ComplexIterator1D       t = target->begin();
        const ComplexIterator1D e = target->end();
        for ( ; t!=e ; ++t )
            *t = x;
    }
}

/** Copy a scalar to all elements. */
void copy(Numeric x,
          ComplexIterator2D target,
          const ComplexIterator2D& end)
{
  for ( ; target!=end ; ++target )
    {
        ComplexIterator1D       t = target->begin();
        const ComplexIterator1D e = target->end();
      for ( ; t!=e ; ++t )
        *t = x;
    }
}


// Functions for ComplexMatrix:
// ---------------------

/** Default constructor. */
ComplexMatrix::ComplexMatrix() :
ComplexMatrixView::ComplexMatrixView()
{
  // Nothing to do here. However, note that the default constructor
  // for ComplexMatrixView has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized. 
}

/** Constructor setting size. This constructor has to set the stride
    in the row range correctly! */
ComplexMatrix::ComplexMatrix(Index r, Index c) :
ComplexMatrixView( new Complex[r*c],
                   ComplexRange(0,r,c),
                   ComplexRange(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
ComplexMatrix::ComplexMatrix(Index r, Index c, Complex fill) :
ComplexMatrixView( new Complex[r*c],
                   ComplexRange(0,r,c),
                   ComplexRange(0,c))
{
    // Here we can access the raw memory directly, for slightly
    // increased efficiency:
    const Complex *stop = mdata+r*c;
    for ( Complex *x=mdata; x<stop; ++x )
        *x = fill;
}

/** Constructor setting size and filling with constant value. */
ComplexMatrix::ComplexMatrix(Index r, Index c, Numeric fill) :
ComplexMatrixView( new Complex[r*c],
              ComplexRange(0,r,c),
            ComplexRange(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
    const Complex *stop = mdata+r*c;
    for ( Complex *x=mdata; x<stop; ++x )
    *x = fill;
}

/** Copy constructor from MatrixView. This automatically sets the size
    and copies the data. */
ComplexMatrix::ComplexMatrix(const ConstComplexMatrixView& m) :
ComplexMatrixView( new Complex[m.nrows()*m.ncols()],
                   ComplexRange( 0, m.nrows(), m.ncols() ),
                   ComplexRange( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Matrix. This automatically sets the size
    and copies the data. */
ComplexMatrix::ComplexMatrix(const ComplexMatrix& m) :
ComplexMatrixView( new Complex[m.nrows()*m.ncols()],
            ComplexRange( 0, m.nrows(), m.ncols() ),
            ComplexRange( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty matrix, then it will have
  // 0 colunns. But this is used to initialize the stride of the row
  // Range! Thus, this method has to be consistent with the behaviour
  // of Range::Range. For now, Range::Range allows also stride 0.
  copy(m.begin(),m.end(),begin());
}

//! Assignment operator from another matrix.
/*! 
  While dimensions of MatrixViews can not be adjusted, dimensions of
  matrices *can* be adjusted. Hence, the behavior of the assignment
  operator is different.

  In this case the size of the target is automatically adjusted. This
  is important, so that structures containing matrices are copied
  correctly. 
  
  This is a deviation from the old ARTS paradigm that sizes must match
  exactly before copying!

  Note: It is sufficient to have only this one version of the
  assignment (Matrix = Matrix). It implicitly covers the cases
  Matrix=MatrixView, etc, because there is a default constructor for
  Matrix from MatrixView. (See C++ Primer Plus, page 571ff.)

  \param m The other matrix to assign to this one.

  \return This matrix, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
ComplexMatrix& ComplexMatrix::operator=(ComplexMatrix m)
{
  swap(*this, m);
  return *this;
}

/** Assignment operator from scalar. Assignment operators also seem to
    be not inherited. */
ComplexMatrix& ComplexMatrix::operator=(Complex x)
{
  copy( x, begin(), end() );
  return *this;
}

//! Assignment from a vector. 
/*! 
  This copies the data from a VectorView to this MatrixView.

  The dimension is adjusted automatically.

  \param v The vector to assign to this matrix.

  \return This matrix, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
ComplexMatrix& ComplexMatrix::operator=(const ConstComplexVectorView& v)
{
  resize( v.nelem(), 1 ); 
  ConstComplexMatrixView dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Matrix is not
    initialized, so it will contain random values.*/
void ComplexMatrix::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );

  if ( mrr.mextent!=r || mcr.mextent!=c )
    {
      delete[] mdata;
      mdata = new Complex[r*c];

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
}


/** Swaps two objects. */
void swap(ComplexMatrix& m1, ComplexMatrix& m2)
{
  std::swap(m1.mrr, m2.mrr);
  std::swap(m1.mcr, m2.mcr);
  std::swap(m1.mdata, m2.mdata);
}


/** Destructor for Matrix. This is important, since Matrix uses new to
    allocate storage. */
ComplexMatrix::~ComplexMatrix()
{
//   cout << "Destroying a Matrix:\n"
//        << *this << "\n........................................\n";
  delete[] mdata;
}


/** Const version of transpose. */
ConstComplexMatrixView transpose(ConstComplexMatrixView m)
{
    return ConstComplexMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
    transpose. The original is not changed! */
ComplexMatrixView transpose(ComplexMatrixView m)
{
    return ComplexMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
 transpose. The original is not changed! */
ComplexMatrixView transpose(ComplexVector v)
{
    return transpose((ComplexMatrixView)v);
}



/** Assignment operator from Array<Complex>. This copies the data from
    an Array<Complex> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Complex> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy. */
ComplexVectorView& ComplexVectorView::operator=(const Array<Complex>& v)
{
  //  cout << "Assigning VectorView from Array<Numeric>.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.nelem());

  // Iterators for Array:
  Array<Complex>::const_iterator i=v.begin();
  const Array<Complex>::const_iterator e=v.end();
  // Iterator for Vector:
  ComplexIterator1D target = begin();

  for ( ; i!=e ; ++i,++target )
    *target = *i;

  return *this;
}

//Functions operating on the complex vectors

/** Scalar product. The two vectors may be identical. */
Complex operator*(const ConstComplexVectorView& a, const ConstComplexVectorView& b)
{
    // Check dimensions:
    assert( a.nelem() == b.nelem() );
    
    const ConstComplexIterator1D ae = a.end();
    ConstComplexIterator1D       ai = a.begin();
    ConstComplexIterator1D       bi = b.begin();
    
    Complex res = 0;
    for ( ; ai!=ae ; ++ai, ++bi )
        res += (*ai) * (*bi);
    
    return res;
}

//! Matrix-Vector Multiplication
/*!
 * 
 * Computes the Matrix-Vector product y = M * x, for a m times n matrix M, a
 * length-m vector y and a length-n vector x.
 * 
 * If you need this to be faster, you have to make an implementation that interacts
 * with cgemv_ in blas.
 * 
 * No memory is allocated for the computation and the matrix and vector views
 * may not overlap.
 * 
 * \param[out] y The length-m ComplexVectorView where the result is stored.
 * \param[in] M Reference to the m-times-n Const{Complex,}MatrixView holding the matrix M.
 * \param[in] x Reference to the length-n Const{Complex,}VectorView holding the vector x.
 */
void mult( ComplexVectorView y,
           const ConstComplexMatrixView& M,
           const ConstComplexVectorView& x )
{
    assert( y.mrange.get_extent() == M.mrr.get_extent() );
    assert( M.mcr.get_extent() == x.mrange.get_extent() );
    assert( (M.mcr.get_extent() != 0) && (M.mrr.get_extent() != 0));
    
    if ((M.mcr.get_stride() == 1) || (M.mrr.get_stride() == 1))
    {
        char trans;
        int m,n;
        std::complex<double> zero = 0.0;
        std::complex<double> one = 1.0;
        int LDA, incx, incy;
        
        if (M.mcr.get_stride() != 1)
        {
            trans = 'n';
            m = (int) M.mrr.get_extent();
            n = (int) M.mcr.get_extent();
            LDA = (int) M.mcr.get_stride();
        }
        else
        {
            trans = 't';
            m = (int) M.mcr.get_extent();
            n = (int) M.mrr.get_extent();
            LDA = (int) M.mrr.get_stride();
            if (M.mrr.get_stride() == 1)
                LDA = m;
        }
        
        incx = (int) x.mrange.get_stride();
        incy = (int) y.mrange.get_stride();
        
        std::complex<double> *mstart = M.mdata + M.mcr.get_start() + M.mrr.get_start();
        std::complex<double> *ystart = y.mdata + y.mrange.get_start();
        std::complex<double> *xstart = x.mdata + x.mrange.get_start();
        
        zgemv_( &trans, &m, &n, &one, mstart, &LDA,
                xstart, &incx, &zero, ystart, &incy );
        
    }
    else
    {
        mult_general( y, M, x );
    }
    
}

/** ComplexMatrix ComplexVector multiplication. y = M*x. Note that the order is different
 *   from MTL, output comes first! Dimensions of y, M, and x must
 *   match. No memory reallocation takes place, only the data is
 *   copied. Using this function on overlapping ComplexMatrix and ComplexVectorViews belonging
 *   to the same Matrix will lead to unpredictable results.
 * 
 *   The implementation here is different from the other multiplication
 *   routines. It does not use iterators but a more drastic approach to gain
 *   maximum performance.  */
void mult_general( ComplexVectorView y,
                   const ConstComplexMatrixView& M,
                   const ConstComplexVectorView& x )
{
    
    // Check dimensions:
    assert( y.mrange.mextent == M.mrr.mextent );
    assert( M.mcr.mextent == x.mrange.mextent );
    assert( M.mcr.mextent != 0 && M.mrr.mextent != 0);
    
    // Let's first find the pointers to the starting positions
    Complex *mdata   = M.mdata + M.mcr.mstart + M.mrr.mstart;
    Complex *xdata   = x.mdata + x.mrange.mstart;
    Complex *yelem   = y.mdata + y.mrange.mstart;
    
    Index i = M.mrr.mextent;
    while (i--)
    {
        Complex *melem = mdata;
        Complex *xelem = xdata; // Reset xelem to first element of source vector
        
        // Multiply first element of matrix row with first element of source
        // vector. This is done outside the loop to avoid initialization of the
        // target vector's element with zero (saves one assignment)
        *yelem = *melem * *xelem;
        
        Index j = M.mcr.mextent;   // --j (instead of j-- like in the outer loop)
        while (--j)                // is important here because we only want
        {                        // mextent-1 cycles
            melem += M.mcr.mstride;
            xelem += x.mrange.mstride;
            *yelem += *melem * *xelem;
        }
        
        mdata += M.mrr.mstride;     // Jump to next matrix row
        yelem += y.mrange.mstride;  // Jump to next element in target vector
    }
}

//! Matrix-Matrix Multiplication
/*!
 * Performs the matrix multiplication A = B * C. The dimensions must match, i.e.
 * A must be a m times n matrix, B a m times k matrix and C a k times c matrix.
 * No memory reallocation takes place, only the data is copied. Using this function
 * on overlapping ComplexMatrixViews belonging to the same ComplexMatrix will lead to unpredictable
 * results. In particular, this means that A and B must not be the same matrix!
 * 
 * If the memory layout allows it, the multiplication is performed using BLAS'
 * _zgemm, which leads to a significant speed up of the operation. To be compatible
 * with BLAS the matrix views A, B and C must satisfy:
 * 
 * - A must have column or row stride 1
 * - B must have column or row stride 1
 * - C must have column stride 1
 * 
 * That means that A and B can be ConstComplexMatrixView objects corresponding to
 * transposed/non-transposed submatrices of a matrix, that are continuous along their
 * first/second dimension. C must correspond to a non-transposed submatrix of a
 * matrix, that is continuous along its second dimension.
 * 
 * Note that additional functions for non-complex complex multiplications must
 * be added as separate routines.
 * 
 * \param[in,out] A The matrix A, that will hold the result of the multiplication.
 * \param[in] B The matrix B
 * \param[in] C The matrix C
 */
void mult( ComplexMatrixView A,
           const ConstComplexMatrixView& B,
           const ConstComplexMatrixView& C )
{
    
    // Check dimensions:
    assert( A.nrows() == B.nrows() );
    assert( A.ncols() == C.ncols() );
    assert( B.ncols() == C.nrows() );
    
    // Catch trivial case if one of the matrices is empty.
    if ( (B.nrows() == 0) || (B.ncols() == 0) || (C.ncols() == 0) )
        return;
    
    // Matrices B and C must be continuous in at least on dimension,  C
    // must be continuous along the second dimension.
    if ( ((B.mrr.get_stride() == 1) || (B.mcr.get_stride() == 1)) &&
        ((C.mrr.get_stride() == 1) || (C.mcr.get_stride() == 1)) &&
        (A.mcr.get_stride() == 1) )
    {
        // BLAS uses column-major order while arts uses row-major order.
        // Hence instead of C = A * B we compute C^T = A^T * B^T!
        
        int k, m, n;
        
        k = (int) B.ncols();
        m = (int) C.ncols();
        n = (int) B.nrows();
        
        // Note also the clash in nomenclature: BLAS uses C = A * B while
        // arts uses A = B * C. Taking into accout this and the difference in
        // memory layouts, we need to map the MatrixViews A, B and C to the BLAS
        // arguments as follows:
        // A (arts) -> C (BLAS)
        // B (arts) -> B (BLAS)
        // C (arts) -> A (BLAS)
        
        // Char indicating whether A (BLAS) or B (BLAS) should be transposed.
        char transa, transb;
        // Sizes of the matrices along the direction in which they are
        // traversed.
        int lda, ldb, ldc;
        
        // Check if C (arts) is transposed.
        if (C.mrr.get_stride() == 1)
        {
            transa = 'T';
            lda = (int) C.mcr.get_stride();
        } else {
            transa = 'N';
            lda = (int) C.mrr.get_stride();
        }
        
        // Check if B (arts) is transposed.
        if (B.mrr.get_stride() == 1)
        {
            transb = 'T';
            ldb = (int) B.mcr.get_stride();
        } else {
            transb = 'N';
            ldb = (int) B.mrr.get_stride();
        }
        
        // In case B (arts) has only one column, column and row stride are 1.
        // We therefore need to set ldb to k, since dgemm_ requires lda to be at
        // least k / m if A is non-transposed / transposed.
        if ( (B.mcr.get_stride() == 1) && (B.mrr.get_stride() == 1) )
        {
            transb = 'N';
            ldb = k;
        }
        
        // The same holds for C (arts).
        if ( (C.mcr.get_stride() == 1) && (C.mrr.get_stride() == 1) )
        {
            transa = 'N';
            lda = m;
        }
        
        ldc = (int) A.mrr.get_stride();
        // The same holds for A (arts).
        if ( (A.mcr.get_stride() == 1) && (A.mrr.get_stride() == 1) )
        {
            ldc = m;
        }
        std::complex<double> alpha = 1.0, beta = 0.0;
        
        zgemm_( & transa,
                & transb,
                & m,
                & n,
                & k,
                & alpha,
                C.mdata + C.mrr.get_start() + C.mcr.get_start(),
                & lda,
                B.mdata + B.mrr.get_start() + B.mcr.get_start(),
                & ldb,
                & beta,
                A.mdata + A.mrr.get_start() + A.mcr.get_start(),
                & ldc );
        
    } else {
        mult_general( A, B, C );
    }
}

//! General matrix multiplication.
/*!
 * This is the fallback matrix multiplication which works for all
 * ConstMatrixView objects.
 * 
 * \param[in,out] A The matrix A, that will hold the result of the multiplication.
 * \param[in] B The matrix B
 * \param[in] C The matrix C
 */
void mult_general( ComplexMatrixView A,
                   const ConstComplexMatrixView& B,
                   const ConstComplexMatrixView& C )
{
    // Check dimensions:
    assert( A.nrows() == B.nrows() );
    assert( A.ncols() == C.ncols() );
    assert( B.ncols() == C.nrows() );
    
    // Let's get the transpose of C, so that we can use 2D iterators to
    // access the columns (= rows of the transpose).
    ConstComplexMatrixView CT = transpose(C);
    
    const ComplexIterator2D ae = A.end();
    ComplexIterator2D       ai = A.begin();
    ConstComplexIterator2D  bi = B.begin();
    
    // This walks through the rows of A and B:
    for ( ; ai!=ae ; ++ai, ++bi )
    {
        const ComplexIterator1D ace = ai->end();
        ComplexIterator1D       aci = ai->begin();
        ConstComplexIterator2D  cti = CT.begin();
        
        // This walks through the columns of A with a 1D iterator, and
        // at the same time through the rows of CT, which are the columns of
        // C, with a 2D iterator:
        for ( ; aci!=ace ; ++aci, ++cti )
        {
            // The operator * is used to compute the scalar product
            // between rows of B and rows of C.transpose().
            *aci = (*bi) * (*cti);
        }
    }
}

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

/** Helper function to access matrix elements.

    Because of function inlining the operator() is not
    accessible from the debuggger. This function helps to access
    Matrix elements from within the debugger.

    \param mv MatrixView
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Complex debug_matrixview_get_elem (ComplexMatrixView& mv, Index r, Index c)
{
  return mv(r, c);
}

#endif
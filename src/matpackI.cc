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
   \file   matpackI.cc

   \author Stefan Buehler
   \date   2001-09-15
*/

#include "blas.h"
#include <cstring>
#include <cmath>
#include "matpackI.h"
#include "exceptions.h"

using std::cout;
using std::endl;
using std::runtime_error;
using std::setw;

// Define the global joker object:
extern const Joker joker = Joker();

static const Numeric RAD2DEG = 57.295779513082323;

// Functions for Range:
// --------------------

/* Explicit constructor. 

  \param Start must be >= 0.

  \param Extent also. Although internally negative extent means "to the end",
  this can not be created this way, only with the joker. Zero
  extent is allowed, though, which corresponds to an empty range.

  \param Stride can be anything. It can be omitted, in which case the
  default value is 1. */
Range::Range(Index start, Index extent, Index stride) :
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
    this means "to the end", or "to the beginning". */
Range::Range(Index start, Joker, Index stride) :
  mstart(start), mextent(-1), mstride(stride)
{
  // Start must be >= 0:
  assert( 0<=mstart );
}

/** Constructor with just a joker. This means, take everything. You
    can still optionally give a stride, though. This constructor is
    just shorter notation for Range(0,joker) */
Range::Range(Joker, Index stride) :
  mstart(0), mextent(-1), mstride(stride)
{
  // Nothing to do here.
}

/** Constructor which converts a range with joker to an explicit
    range.

    \param max_size The maximum allowed size of the vector. 
    \param r The new range, with joker. */
Range::Range(Index max_size, const Range& r) :
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
    may contain -1 for the stride, which acts as a joker.

    \param p Previous range.
    \param n New range. */
Range::Range(const Range& p, const Range& n) :
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
  // and <= prev_fin, except for Joker:
  assert( mstart<=prev_fin || mextent == -1 );

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

// Functions for ConstVectorView:
// ------------------------------

//! Returns true if variable size is zero.
bool ConstVectorView::empty() const
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
Index ConstVectorView::nelem() const
{
  return mrange.mextent;
}

/** The sum of all elements of a Vector. */
Numeric ConstVectorView::sum() const
{
  Numeric s=0;
  ConstIterator1D i = begin();
  const ConstIterator1D e = end();

  for ( ; i!=e; ++i )
    s += *i;

  return s;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Vector. This allows
    correct recursive behavior.  */
ConstVectorView ConstVectorView::operator[](const Range& r) const
{
  return ConstVectorView(mdata, mrange, r);
}

/** Return const iterator to first element. */
ConstIterator1D ConstVectorView::begin() const
{
  return ConstIterator1D(mdata+mrange.mstart, mrange.mstride);
}

/** Return const iterator behind last element. */
ConstIterator1D ConstVectorView::end() const
{
  return ConstIterator1D( mdata +
                          mrange.mstart +
                          (mrange.mextent)*mrange.mstride,
                          mrange.mstride );
}

/** Conversion to const 1 column matrix. */
ConstVectorView::operator ConstMatrixView() const
{
    // The old version (before 2013-01-18) of this was:
    //    return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
    // Bus this was a bug! The problem is that the matrix index operator adds
    // the mstart from both row and columm range object to mdata

    return ConstMatrixView(mdata,mrange,Range(0,1));
}

/** A special constructor, which allows to make a ConstVectorView from
    a scalar.

    This one is a bit tricky: We have to cast away the arguments const
    qualifier, because mdata is not const. This should be safe, since
    there are no non-const methods for ConstVectorView.
*/
ConstVectorView::ConstVectorView(const Numeric& a) :
  mrange(0,1), mdata(&const_cast<Numeric&>(a))
{
  // Nothing to do here.
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Vector. */
ConstVectorView::ConstVectorView() :
  mrange(0,0), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
ConstVectorView::ConstVectorView( Numeric *data,
                                  const Range& range) :
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
ConstVectorView::ConstVectorView( Numeric *data,
                                  const Range& p,
                                  const Range& n) :
  mrange(p,n),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the vector. The iterators know which part of the vector
    is `active', and also the stride. */
std::ostream& operator<<(std::ostream& os, const ConstVectorView& v)
{
  ConstIterator1D i=v.begin();
  const ConstIterator1D end=v.end();

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


// Functions for VectorView:
// ------------------------

/** Bail out immediately if somebody tries to create a VectorView from
    a const Vector. */
VectorView::VectorView (const Vector&)
{
  throw runtime_error("Creating a VectorView from a const Vector is not allowed.");
  // This is not really a runtime error, but I don't want to start
  // producing direct output from inside matpack. And just exiting is
  // not so nice.
  // If you see this error, there is a bug in the code, not in the
  // ARTS input.
}

/** Create VectorView from a Vector. */
VectorView::VectorView (Vector& v)
{
  mdata = v.mdata;
  mrange = v.mrange;
}

/** Const index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  Has to be redifined here, because the
    one from ConstVectorView is hidden. */
ConstVectorView VectorView::operator[](const Range& r) const
{
  return ConstVectorView::operator[](r);
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  */
VectorView VectorView::operator[](const Range& r)
{
  return VectorView(mdata, mrange, r);
}

/** Return const iterator to first element. Has to be redefined here,
    since it is hiden by the non-const operator of the derived
    class.*/
ConstIterator1D VectorView::begin() const
{
  return ConstVectorView::begin();
}

/** Return const iterator behind last element. Has to be redefined
    here, since it is hiden by the non-const operator of the derived
    class.*/
ConstIterator1D VectorView::end() const
{
  return ConstVectorView::end();
}

/** Return iterator to first element. */
Iterator1D VectorView::begin()
{
  return Iterator1D(mdata+mrange.mstart, mrange.mstride);
}

/** Return iterator behind last element. */
Iterator1D VectorView::end()
{
  return Iterator1D( mdata +
                     mrange.mstart +
                     (mrange.mextent)*mrange.mstride,
                     mrange.mstride );
}

/** Assignment operator. This copies the data from another VectorView
    to this VectorView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this VectorView by
    setting its range. */
VectorView& VectorView::operator=(const ConstVectorView& v)
{
  //  cout << "Assigning VectorView from ConstVectorView.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);

  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assignment from VectorView to VectorView. This is a tricky
    one. The problem is that since VectorView is derived from
    ConstVectorView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
VectorView& VectorView::operator=(const VectorView& v)
{
  //  cout << "Assigning VectorView from VectorView.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);

  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assignment from Vector. This is important to avoid a bug when
    assigning a Vector to a VectorView. */
VectorView& VectorView::operator=(const Vector& v)
{
  //  cout << "Assigning VectorView from Vector.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);

  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assigning a scalar to a VectorView will set all elements to this
    value. */
VectorView& VectorView::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
VectorView VectorView::operator*=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i *= x;
  return *this;
}

/** Division by scalar. */
VectorView VectorView::operator/=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i /= x;
  return *this;
}

/** Addition of scalar. */
VectorView VectorView::operator+=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i += x;
  return *this;
}

/** Subtraction of scalar. */
VectorView VectorView::operator-=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i -= x;
  return *this;
}

/** Element-vise multiplication by another vector. */
VectorView VectorView::operator*=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i *= *s;
  return *this;
}

/** Element-vise division by another vector. */
VectorView VectorView::operator/=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i /= *s;
  return *this;
}

/** Element-vise addition of another vector. */
VectorView VectorView::operator+=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i += *s;
  return *this;
}

/** Element-vise subtraction of another vector. */
VectorView VectorView::operator-=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i -= *s;
  return *this;
}

/** Conversion to 1 column matrix. */
VectorView::operator MatrixView()
{
    // The old version (before 2013-01-18) of this was:
    //    return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
    // Bus this was a bug! The problem is that the matrix index operator adds
    // the mstart from both row and columm range object to mdata

    return MatrixView(mdata,mrange,Range(0,1));
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  VectorView is not pointing to the beginning of a Vector or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric * VectorView::get_c_array() const
{
    if (mrange.mstart != 0 || mrange.mstride != 1)
        throw std::runtime_error("A VectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");

    return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  VectorView is not pointing to the beginning of a Vector or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric * VectorView::get_c_array()
{
    if (mrange.mstart != 0 || mrange.mstride != 1)
        throw std::runtime_error("A VectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** A special constructor, which allows to make a VectorView from
    a scalar.
*/
VectorView::VectorView(Numeric& a) :
  ConstVectorView(a)
{
  // Nothing to do here.
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Vector. */
VectorView::VectorView() :
  ConstVectorView()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
VectorView::VectorView(Numeric *data,
                            const Range& range) :
  ConstVectorView(data,range)
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
VectorView::VectorView(Numeric *data,
                              const Range& p,
                              const Range& n) :
  ConstVectorView(data,p,n)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors. */
void copy(ConstIterator1D origin,
                 const ConstIterator1D& end,
                 Iterator1D target)
{
  if (origin.mstride == 1 && target.mstride == 1)
    memcpy((void *)target.mx,
           (void *)origin.mx,
           sizeof(Numeric) * (end.mx - origin.mx));
  else
    for ( ; origin!=end ; ++origin,++target )
      *target = *origin;
}

/** Copy a scalar to all elements. */
void copy(Numeric x,
          Iterator1D target,
          const Iterator1D& end)
{
  for ( ; target!=end ; ++target )
    *target = x;
}


// Functions for Vector:
// ---------------------

/** Default constructor. */
Vector::Vector()
{
  // Nothing to do here
}

/** Initialization list constructor. */
Vector::Vector(std::initializer_list<Numeric> init) : VectorView(
        new Numeric[init.size()],
        Range(0, init.size()))
{
  std::copy(init.begin(), init.end(), begin());
}

/** Constructor setting size. */
Vector::Vector(Index n) :
  VectorView( new Numeric[n],
             Range(0,n))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Vector::Vector(Index n, Numeric fill) :
  VectorView( new Numeric[n],
             Range(0,n))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, n, fill);
}

/** Constructor filling with values. 

    Examples:

    Vector v(1,5,1);  // 1, 2, 3, 4, 5
    Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
    Vector v(5,5,-1); // 5, 4, 3, 2, 1
*/
Vector::Vector(Numeric start, Index extent, Numeric stride) :
  VectorView( new Numeric[extent],
             Range(0,extent))
{
  // Fill with values:
  Numeric x = start;
  Iterator1D       i=begin();
  const Iterator1D e=end();
  for ( ; i!=e; ++i )
    {
      *i = x;
      x += stride;
    }
}

/** Copy constructor from VectorView. This automatically sets the size
    and copies the data. The vector created will have start zero and
    stride 1, independent on how these parameters are set for the
    original. So, what is copied is the data, not the shape
    of the selection. */
Vector::Vector(const ConstVectorView& v) :
  VectorView( new Numeric[v.nelem()],
              Range(0,v.nelem()))
{
  copy(v.begin(),v.end(),begin());
}

/** Copy constructor from Vector. This is important to override the
    automatically generated shallow constructor. We want deep copies!  */
Vector::Vector(const Vector& v) :
  VectorView( new Numeric[v.nelem()],
              Range(0,v.nelem()))
{
  std::memcpy(mdata, v.mdata, nelem()*sizeof(Numeric));
}

/** Converting constructor from std::vector<Numeric>. */
Vector::Vector(const std::vector<Numeric>& v) :
  VectorView( new Numeric[v.size()],
              Range(0,v.size()))
{
    std::vector<Numeric>::const_iterator vec_it_end = v.end();
    Iterator1D this_it = this->begin();
    for (std::vector<Numeric>::const_iterator vec_it = v.begin();
         vec_it != vec_it_end;
         ++vec_it, ++this_it)
        *this_it = *vec_it;
}

//! Assignment from an initializatoin list.
Vector& Vector::operator=(std::initializer_list<Numeric> v)
{
  resize(v.size());
  std::copy(v.begin(), v.end(), begin());
  return *this;
}

//! Assignment from another Vector.
/*! 
  While dimensions of VectorViews can not be adjusted, dimensions of
  Vectors *can* be adjusted. Hence, the behavior of the assignment
  operator is different.

  In this case the size of the target is automatically adjusted. This
  is important, so that structures containing Vectors are copied
  correctly. 

  This is a deviation from the old ARTS paradigm that sizes must match
  exactly before copying!

  Note: It is sufficient to have only this one version of the
  assignment (Vector = Vector). It implicitly covers the cases
  Vector=VectorView, etc, because there is a default constructor for
  Vector from VectorView. (See C++ Primer Plus, page 571ff.)

  \param v The other vector to copy to this one.

  \return This vector, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
Vector& Vector::operator=(const Vector& v)
{
  if (this != &v)
  {
    resize(v.nelem());
    std::memcpy(mdata, v.mdata, nelem()*sizeof(Numeric));
  }
  return *this;
}

//! Move assignment from another Vector.
Vector& Vector::operator=(Vector&& v) noexcept
{
  if (this != &v)
  {
    delete[] mdata;
    mdata = v.mdata;
    mrange = v.mrange;
    v.mrange = Range(0, 0);
    v.mdata = nullptr;
  }
  return *this;
}

//! Assignment operator from Array<Numeric>.
/*!
  This copies the data from a Array<Numeric> to this VectorView. The
  size is adjusted automatically.

  Array<Numeric> can be useful to collect things in, because there
  is a .push_back method for it. Then, after collecting we usually
  have to transfer the content to a Vector. With this assignment
  operator that's easy.

  \param x The array to assign to this.

  \return This vector, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
Vector& Vector::operator=(const Array<Numeric>& x)
{
  resize( x.nelem() );
  VectorView::operator=(x);
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
Vector& Vector::operator=(Numeric x)
{
  std::fill_n(mdata, nelem(), x);
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
    nothing. All data is lost after resizing! The new Vector is not
    initialized, so it will contain random values.  */
void Vector::resize(Index n)
{
  assert( 0<=n );
  if ( mrange.mextent != n )
    {
      delete[] mdata;
      mdata = new Numeric[n];
      mrange.mstart = 0;
      mrange.mextent = n;
      mrange.mstride = 1;
    }
}


/** Swaps two objects. */
void swap(Vector& v1, Vector& v2)
{
  std::swap(v1.mrange, v2.mrange);
  std::swap(v1.mdata, v2.mdata);
}


/** Destructor for Vector. This is important, since Vector uses new to
    allocate storage. */
Vector::~Vector()
{
  delete[] mdata;
}


// Functions for ConstMatrixView:
// ------------------------------

//! Returns true if variable size is zero.
bool ConstMatrixView::empty() const
{
    return (nrows() == 0 || ncols() == 0);
}

/** Returns the number of rows. */
Index ConstMatrixView::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
Index ConstMatrixView::ncols() const
{
  return mcr.mextent;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior.  */
ConstMatrixView ConstMatrixView::operator()(const Range& r,
                                                   const Range& c) const
{
  return ConstMatrixView(mdata, mrr, mcr, r, c);
}

/** Const index operator returning a column as an object of type
    ConstVectorView.

    \param r A range of rows.
    \param c Index of selected column */
ConstVectorView ConstMatrixView::operator()(const Range& r, Index c) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata + mcr.mstart + c*mcr.mstride,
                         mrr, r);
}

/** Const index operator returning a row as an object of type
    ConstVectorView.

    \param r Index of selected row.
    \param c Range of columns */
ConstVectorView ConstMatrixView::operator()(Index r, const Range& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return ConstVectorView(mdata + mrr.mstart + r*mrr.mstride,
                         mcr, c);
}

/** Return const iterator to first row. */
ConstIterator2D ConstMatrixView::begin() const
{
  return ConstIterator2D(ConstVectorView(mdata+mrr.mstart,
                                         mcr),
                         mrr.mstride);
}

/** Return const iterator behind last row. */
ConstIterator2D ConstMatrixView::end() const
{
  return ConstIterator2D( ConstVectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
                                          mcr),
                          mrr.mstride );
}

//! Matrix diagonal as vector.
/*!
  Returns a ConstMatrixView on the diagonal entries of the matrix. For a given
  (n,m) matrix M the diagonal vector v is the vector of length min{n,m} with entries

       v[i] = M(i,i)

  \return The diagonal vector v.
*/
ConstVectorView ConstMatrixView::diagonal() const
{
    Index n = std::min( mrr.mextent, mcr.mextent );
    return ConstVectorView( mdata + mrr.mstart + mcr.mstart,
                            Range( 0, n, mrr.mstride + mcr.mstride ) );
}


/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
ConstMatrixView::ConstMatrixView() :
  mrr(0,0,1), mcr(0,0,1), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
ConstMatrixView::ConstMatrixView( Numeric *data,
                                  const Range& rr,
                                  const Range& cr) :
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
ConstMatrixView::ConstMatrixView( Numeric *data,
                                  const Range& pr, const Range& pc,
                                  const Range& nr, const Range& nc) :
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
std::ostream& operator<<(std::ostream& os, const ConstMatrixView& v)
{
  // Row iterators:
  ConstIterator2D ir=v.begin();
  const ConstIterator2D end_row=v.end();

  if ( ir!=end_row )
    {
      ConstIterator1D ic =  ir->begin();
      const ConstIterator1D end_col = ir->end();

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
      ConstIterator1D ic =  ir->begin();
      const ConstIterator1D end_col = ir->end();

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


// Functions for MatrixView:
// -------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
ConstMatrixView MatrixView::operator()(const Range& r, const Range& c) const
{
  return ConstMatrixView::operator()(r,c);
}

/** Const index operator returning a column as an object of type
    ConstVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r A range of rows.
    \param c Index of selected column */
ConstVectorView MatrixView::operator()(const Range& r, Index c) const
{
  return ConstMatrixView::operator()(r,c);
}

/** Const index operator returning a row as an object of type
    ConstVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r Index of selected row.
    \param c Range of columns */
ConstVectorView MatrixView::operator()(Index r, const Range& c) const
{
  return ConstMatrixView::operator()(r,c);
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Matrix. This allows correct
    recursive behavior.  */
MatrixView MatrixView::operator()(const Range& r, const Range& c)
{
  return MatrixView(mdata, mrr, mcr, r, c);
}

/** Index operator returning a column as an object of type VectorView.

    \param r A range of rows.
    \param c Index of selected column */
VectorView MatrixView::operator()(const Range& r, Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return VectorView(mdata + mcr.mstart + c*mcr.mstride,
                    mrr, r);
}

/** Index operator returning a row as an object of type VectorView.

    \param r Index of selected row.
    \param c Range of columns */
VectorView MatrixView::operator()(Index r, const Range& c)
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r <  mrr.mextent );

  return VectorView(mdata + mrr.mstart + r*mrr.mstride,
                    mcr, c);
}

/** Return const iterator to first row. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
ConstIterator2D MatrixView::begin() const
{
  return ConstMatrixView::begin();
}

/** Return const iterator behind last row. */
ConstIterator2D MatrixView::end() const
{
  return ConstMatrixView::end();
}

/** Return iterator to first row. */
Iterator2D MatrixView::begin()
{
  return Iterator2D(VectorView(mdata+mrr.mstart, mcr),
                    mrr.mstride);
}

/** Return iterator behind last row. */
Iterator2D MatrixView::end()
{
  return Iterator2D( VectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
                                mcr),
                     mrr.mstride );
}

/** Assignment operator. This copies the data from another MatrixView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
MatrixView& MatrixView::operator=(const ConstMatrixView& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from MatrixView to MatrixView. This is a tricky
    one. The problem is that since MatrixView is derived from
    ConstMatrixView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
MatrixView& MatrixView::operator=(const MatrixView& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Matrix. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
MatrixView& MatrixView::operator=(const Matrix& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a vector. This copies the data from a VectorView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
MatrixView& MatrixView::operator=(const ConstVectorView& v)
{
  // Check that sizes are compatible:
  assert( mrr.mextent==v.nelem() );
  assert( mcr.mextent==1         );
  //  dummy = ConstMatrixView(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstMatrixView dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Assigning a scalar to a MatrixView will set all elements to this
    value. */
MatrixView& MatrixView::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
MatrixView& MatrixView::operator*=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c *= x;
    }
  return *this;
}

/** Division by scalar. */
MatrixView& MatrixView::operator/=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c /= x;
    }
  return *this;
}

/** Addition of scalar. */
MatrixView& MatrixView::operator+=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c += x;
    }
  return *this;
}

/** Subtraction of scalar. */
MatrixView& MatrixView::operator-=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c -= x;
    }
  return *this;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric *MatrixView::get_c_array() const
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
Numeric *MatrixView::get_c_array()
{
    if (mrr.mstart != 0 || mrr.mstride != mcr.mextent
        || mcr.mstart != 0 || mcr.mstride != 1)
        throw std::runtime_error("A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");

    return mdata;
}

/** Element-vise multiplication by another Matrix. */
MatrixView& MatrixView::operator*=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin();
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c *= *sc;
    }
  return *this;
}

/** Element-vise division by another Matrix. */
MatrixView& MatrixView::operator/=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin();
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c /= *sc;
    }
  return *this;
}

/** Element-vise addition of another Matrix. */
MatrixView& MatrixView::operator+=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin();
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c += *sc;
    }
  return *this;
}

/** Element-vise subtraction of another Matrix. */
MatrixView& MatrixView::operator-=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin();
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c -= *sc;
    }
  return *this;
}

/** Element-vise multiplication by a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator*=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c *= *sc;
    }
  return *this;
}

/** Element-vise division by a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator/=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c /= *sc;
    }
  return *this;
}

/** Element-vise addition of a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator+=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c += *sc;
    }
  return *this;
}

/** Element-vise subtraction of a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator-=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c -= *sc;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Matrix. */
MatrixView::MatrixView() :
  ConstMatrixView()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
MatrixView::MatrixView(Numeric *data,
                            const Range& rr,
                            const Range& cr) :
  ConstMatrixView(data, rr, cr)
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
MatrixView::MatrixView(Numeric *data,
                            const Range& pr, const Range& pc,
                            const Range& nr, const Range& nc) :
  ConstMatrixView(data,pr,pc,nr,nc)
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
void copy(ConstIterator2D origin,
                 const ConstIterator2D& end,
                 Iterator2D target)
{
  for ( ; origin!=end ; ++origin,++target )
    {
      ConstIterator1D       o = origin->begin();
      const ConstIterator1D e = origin->end();
      Iterator1D            t = target->begin();
      for ( ; o!=e ; ++o,++t )
        *t = *o;
    }
}

/** Copy a scalar to all elements. */
void copy(Numeric x,
                 Iterator2D target,
                 const Iterator2D& end)
{
  for ( ; target!=end ; ++target )
    {
      Iterator1D       t = target->begin();
      const Iterator1D e = target->end();
      for ( ; t!=e ; ++t )
        *t = x;
    }
}


// Functions for Matrix:
// ---------------------

/** Default constructor. */
Matrix::Matrix() :
  MatrixView::MatrixView()
{
  // Nothing to do here. However, note that the default constructor
  // for MatrixView has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized.
}

/** Constructor setting size. This constructor has to set the stride
    in the row range correctly! */
Matrix::Matrix(Index r, Index c) :
  MatrixView( new Numeric[r*c],
             Range(0,r,c),
             Range(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Matrix::Matrix(Index r, Index c, Numeric fill) :
  MatrixView( new Numeric[r*c],
              Range(0,r,c),
              Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, r*c, fill);
}

/** Copy constructor from MatrixView. This automatically sets the size
    and copies the data. */
Matrix::Matrix(const ConstMatrixView& m) :
  MatrixView( new Numeric[m.nrows()*m.ncols()],
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Matrix. This automatically sets the size
    and copies the data. */
Matrix::Matrix(const Matrix& m) :
  MatrixView( new Numeric[m.nrows()*m.ncols()],
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty matrix, then it will have
  // 0 colunns. But this is used to initialize the stride of the row
  // Range! Thus, this method has to be consistent with the behaviour
  // of Range::Range. For now, Range::Range allows also stride 0.
  std::memcpy(mdata, m.mdata, nrows()*ncols()*sizeof(Numeric));
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
Matrix& Matrix::operator=(const Matrix& m)
{
  if (this != &m)
  {
    resize(m.nrows(), m.ncols());
    std::memcpy(mdata, m.mdata, nrows()*ncols()*sizeof(Numeric));
  }
  return *this;
}

//! Move assignment operator from another matrix.
Matrix& Matrix::operator=(Matrix&& m) noexcept
{
  if (this != &m)
  {
    delete[] mdata;
    mdata = m.mdata;
    mrr = m.mrr;
    mcr = m.mcr;
    m.mrr = Range(0, 0);
    m.mcr = Range(0, 0);
    m.mdata = nullptr;
  }
  return *this;
}

/** Assignment operator from scalar. Assignment operators also seem to
    be not inherited. */
Matrix& Matrix::operator=(Numeric x)
{
  std::fill_n(mdata, nrows()*ncols(), x);
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
Matrix& Matrix::operator=(const ConstVectorView& v)
{
  resize( v.nelem(), 1 );
  ConstMatrixView dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Matrix is not
    initialized, so it will contain random values.*/
void Matrix::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );

  if ( mrr.mextent!=r || mcr.mextent!=c )
    {
      delete[] mdata;
      mdata = new Numeric[r*c];

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
}


/** Swaps two objects. */
void swap(Matrix& m1, Matrix& m2)
{
  std::swap(m1.mrr, m2.mrr);
  std::swap(m1.mcr, m2.mcr);
  std::swap(m1.mdata, m2.mdata);
}


/** Destructor for Matrix. This is important, since Matrix uses new to
    allocate storage. */
Matrix::~Matrix()
{
//   cout << "Destroying a Matrix:\n"
//        << *this << "\n........................................\n";
  delete[] mdata;
}


// Some general Matrix Vector functions:

/** Scalar product. The two vectors may be identical. */
Numeric operator*(const ConstVectorView& a, const ConstVectorView& b)
{
  // Check dimensions:
  assert( a.nelem() == b.nelem() );

  const ConstIterator1D ae = a.end();
  ConstIterator1D       ai = a.begin();
  ConstIterator1D       bi = b.begin();

  Numeric res = 0;
  for ( ; ai!=ae ; ++ai, ++bi )
    res += (*ai) * (*bi);

  return res;
}

//! Matrix-Vector Multiplication
/*!

  Computes the Matrix-Vector product y = M * x, for a m times n matrix M, a
  length-m vector y and a length-n vector x.

  The product is computed using the dgemv_ routine from the BLAS library if
  the matrix is contiguous in memory. If this is not the case, the mult_general
  method is used to compute the product.

  No memory is allocated for the computation and the matrix and vector views
  may not overlap.

  \param[out] y The length-m VectorView where the result is stored.
  \param[in] M Reference to the m-times-n ConstMatrixView holding the matrix M.
  \param[in] x Reference to the length-n ConstVectorView holding the vector x.
*/
void mult( VectorView y,
           const ConstMatrixView& M,
           const ConstVectorView& x )
{

    assert( y.mrange.get_extent() == M.mrr.get_extent() );
    assert( M.mcr.get_extent() == x.mrange.get_extent() );
    assert( (M.mcr.get_extent() != 0) && (M.mrr.get_extent() != 0));

    if ((M.mcr.get_stride() == 1) || (M.mrr.get_stride() == 1))
    {
        char trans;
        int m,n;
        double zero = 0.0;
        double one = 1.0;
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

        double *mstart = M.mdata + M.mcr.get_start() + M.mrr.get_start();
        double *ystart = y.mdata + y.mrange.get_start();
        double *xstart = x.mdata + x.mrange.get_start();

        dgemv_( &trans, &m, &n, &one, mstart, &LDA,
                xstart, &incx, &zero, ystart, &incy );

    }
    else
    {
        mult_general( y, M, x );
    }

}

/** Matrix Vector multiplication. y = M*x. Note that the order is different
    from MTL, output comes first! Dimensions of y, M, and x must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping Matrix and VectorViews belonging
    to the same Matrix will lead to unpredictable results.

    The implementation here is different from the other multiplication
    routines. It does not use iterators but a more drastic approach to gain
    maximum performance.  */
void mult_general( VectorView y,
                   const ConstMatrixView& M,
                   const ConstVectorView& x )
{

  // Check dimensions:
  assert( y.mrange.mextent == M.mrr.mextent );
  assert( M.mcr.mextent == x.mrange.mextent );
  assert( M.mcr.mextent != 0 && M.mrr.mextent != 0);

  // Let's first find the pointers to the starting positions
  Numeric *mdata   = M.mdata + M.mcr.mstart + M.mrr.mstart;
  Numeric *xdata   = x.mdata + x.mrange.mstart;
  Numeric *yelem   = y.mdata + y.mrange.mstart;

  Index i = M.mrr.mextent;
  while (i--)
    {
      Numeric *melem = mdata;
      Numeric *xelem = xdata; // Reset xelem to first element of source vector

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
  Performs the matrix multiplication A = B * C. The dimensions must match, i.e.
  A must be a m times n matrix, B a m times k matrix and C a k times c matrix.
  No memory reallocation takes place, only the data is copied. Using this function
  on overlapping MatrixViews belonging to the same Matrix will lead to unpredictable
  results. In particular, this means that A and B must not be the same matrix!

  If the memory layout allows it, the multiplication is performed using BLAS'
  _dgemm, which leads to a significant speed up of the operation. To be compatible
  with BLAS the matrix views A, B and C must satisfy:

  - A must have column or row stride 1
  - B must have column or row stride 1
  - C must have column stride 1

  That means that A and B can be ConstMatrixView objects corresponding to
  transposed/non-transposed submatrices of a matrix, that are continuous along their
  first/second dimension. C must correspond to a non-transposed submatrix of a
  matrix, that is continuous along its second dimension.

  \param[in,out] A The matrix A, that will hold the result of the multiplication.
  \param[in] B The matrix B
  \param[in] C The matrix C
*/
void mult( MatrixView A,
           const ConstMatrixView& B,
           const ConstMatrixView& C )
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
        double alpha = 1.0, beta = 0.0;

        dgemm_( & transa,
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
  This is the fallback matrix multiplication which works for all
  ConstMatrixView objects.

  \param[in,out] A The matrix A, that will hold the result of the multiplication.
  \param[in] B The matrix B
  \param[in] C The matrix C
*/
void mult_general( MatrixView A,
                   const ConstMatrixView& B,
                   const ConstMatrixView& C )
{
  // Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Let's get the transpose of C, so that we can use 2D iterators to
  // access the columns (= rows of the transpose).
  ConstMatrixView CT = transpose(C);

  const Iterator2D ae = A.end();
  Iterator2D       ai = A.begin();
  ConstIterator2D  bi = B.begin();

  // This walks through the rows of A and B:
  for ( ; ai!=ae ; ++ai, ++bi )
    {
      const Iterator1D ace = ai->end();
      Iterator1D       aci = ai->begin();
      ConstIterator2D  cti = CT.begin();

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

//! cross3
/*!
    Calculates the cross product between two vectors of length 3.

    c = a x b, for 3D vectors. The vector c must have length 3 and can not be
    the same variable as a or b.

    param    c   Out: The cross product vector
    \param   a   In: A vector of length 3.
    \param   b   In: A vector of length 3.

    \author Patrick Eriksson
    \date   2012-02-12
*/
void cross3(VectorView c, const ConstVectorView& a, const ConstVectorView& b)
{
  assert( a.nelem() == 3 );
  assert( b.nelem() == 3 );
  assert( c.nelem() == 3 );

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/*!
    Returns numeric angle between two vectors in degrees.

    \param  a   In:  A vector of length N.
    \param  b   In:  A vector of length N.

    \author Richard Larsson
    \date   2012-07-10
*/
Numeric vector_angle(ConstVectorView a, ConstVectorView b)
{
    assert( a.nelem() == b.nelem() );
    Numeric arg = (a*b) / sqrt(a*a) / sqrt(b*b);

    // Due to numerical inaccuracy, arg might be slightly larger than 1.
    // We catch those cases to avoid spurious returns of NaNs
    return fabs(arg) > 1. ? 0. : acos(arg) * RAD2DEG;
}


/*!
    Calculates the projection of two vectors of equal length.

    c = proj_a(b). Projecting b on a. The vector c must have the same length
    but can not be the same variable as a or b.

    \param   c   Out: The projection of b on a.
    \param   a   In:  A vector of length N.
    \param   b   In:  A vector of length N.

    \author Richard Larsson
    \date   2012-07-10
*/
void proj(Vector& c, ConstVectorView a, ConstVectorView b)
{
    assert( a.nelem() == b.nelem() );
    assert( a.nelem() == c.nelem() );

    const Numeric C = (a*b) / (a*a);
        c  = a;
        c *= C;
};


/** Const version of transpose. */
ConstMatrixView transpose(ConstMatrixView m)
{
  return ConstMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
    transpose. The original is not changed! */
MatrixView transpose(MatrixView m)
{
  return MatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
 transpose. The original is not changed! */
MatrixView transpose(Vector v)
{
  return transpose((MatrixView)v);
}

/** A generic transform function for vectors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    Although the matrix version of this can also be used for vectors,
    thanks to the automatic interpretation of a vector as a one column
    matrix, this one is slightly more efficient. However, the
    difference is very small (only a few percent). 

    The two views may be the same one, in which case the
    conversion happens in place. 

    \param   y Output:   The results of the function acting on each element of x.
    \param    my_func A function (e.g., sqrt).
    \param    x   A vector. */
void transform( VectorView y,
                       double (&my_func)(double),
                       ConstVectorView x )
{
  // Check dimensions:
  assert( y.nelem()==x.nelem() );

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();
  Iterator1D            yi = y.begin();
  for ( ; xi!=xe; ++xi, ++yi )
    *yi = my_func(*xi);
}

/** A generic transform function for matrices, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    This function can also be used for Vectors, because there is a
    conversion to MatrixView.

    The two Matrix views may be the same one, in which case the
    conversion happens in place. 

   \param   y Output:   The results of the function acting on each element of x.
   \param    my_func A function (e.g., sqrt).
   \param    x   A matrix. */
void transform( MatrixView y,
                       double (&my_func)(double),
                       ConstMatrixView x )
{
  // Check dimensions:
  assert( y.nrows()==x.nrows() );
  assert( y.ncols()==x.ncols() );

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();
  Iterator2D             ry = y.begin();
  for ( ; rx!=rxe; ++rx, ++ry )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();
      Iterator1D             cy = ry->begin();
      for ( ; cx!=cxe; ++cx, ++cy )
        *cy = my_func(*cx);
    }
}

/** Max function, vector version. */
Numeric max(const ConstVectorView& x)
{
  // Initial value for max:
  Numeric max = x[0];

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();

  for ( ; xi!=xe ; ++xi )
    {
      if ( *xi > max )
        max = *xi;
    }

  return max;
}

/** Max function, matrix version. */
Numeric max(const ConstMatrixView& x)
{
  // Initial value for max:
  Numeric max = x(0,0);

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();

  for ( ; rx!=rxe ; ++rx )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();

      for ( ; cx!=cxe ; ++cx )
        if ( *cx > max )
          max = *cx;
    }

  return max;
}

/** Min function, vector version. */
Numeric min(const ConstVectorView& x)
{
  // Initial value for min:
  Numeric min = x[0];

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();

  for ( ; xi!=xe ; ++xi )
    {
      if ( *xi < min )
        min = *xi;
    }

  return min;
}

/** Min function, matrix version. */
Numeric min(const ConstMatrixView& x)
{
  // Initial value for min:
  Numeric min = x(0,0);

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();

  for ( ; rx!=rxe ; ++rx )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();

      for ( ; cx!=cxe ; ++cx )
        if ( *cx < min )
          min = *cx;
    }

  return min;
}

/** Mean function, vector version. */
Numeric mean(const ConstVectorView& x)
{
  // Initial value for mean:
  Numeric mean = 0;

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();

  for ( ; xi!=xe ; ++xi ) mean += *xi;

  mean /= (Numeric)x.nelem();

  return mean;
}

/** Mean function, matrix version. */
Numeric mean(const ConstMatrixView& x)
{
  // Initial value for mean:
  Numeric mean = 0;

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();

  for ( ; rx!=rxe ; ++rx )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();

      for ( ; cx!=cxe ; ++cx ) mean += *cx;
    }

  mean /= (Numeric)(x.nrows()*x.ncols());

  return mean;
}

/** Assignment operator from Array<Numeric>. This copies the data from
    an Array<Numeric> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Numeric> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy. */
VectorView& VectorView::operator=(const Array<Numeric>& v)
{
  //  cout << "Assigning VectorView from Array<Numeric>.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.nelem());

  // Iterators for Array:
  Array<Numeric>::const_iterator i=v.begin();
  const Array<Numeric>::const_iterator e=v.end();
  // Iterator for Vector:
  Iterator1D target = begin();

  for ( ; i!=e ; ++i,++target )
    *target = *i;

  return *this;
}

// Const

// Converts constant matrix to constant eigen map
ConstMatrixViewMap MapToEigen(const ConstMatrixView& A){
  return ConstMatrixViewMap(A.mdata+A.mrr.get_start()+A.mcr.get_start(),
   A.nrows(), A.ncols(), StrideType(A.mrr.get_stride(), A.mcr.get_stride())); }

// Converts constant vector to constant eigen row-view
ConstMatrixViewMap MapToEigen(const ConstVectorView& A){
  return ConstMatrixViewMap(A.mdata+A.mrange.get_start(),
   A.nelem(), 1, StrideType(A.mrange.get_stride(), 1)); }

// Converts constant vector to constant eigen row-view
ConstMatrixViewMap MapToEigenRow(const ConstVectorView& A){return MapToEigen(A);}

// Converts constant vector to constant eigen column-view
ConstMatrixViewMap MapToEigenCol(const ConstVectorView& A){
  return ConstMatrixViewMap(A.mdata+A.mrange.get_start(),
   1, A.nelem(), StrideType(1, A.mrange.get_stride())); }

// Non- const

// Converts matrix to eigen map
MatrixViewMap MapToEigen(MatrixView& A){
  return MatrixViewMap(A.mdata+A.mrr.get_start()+A.mcr.get_start(),
   A.nrows(), A.ncols(), StrideType(A.mrr.get_stride(), A.mcr.get_stride())); }

// Converts vector to eigen map row-view
MatrixViewMap MapToEigen(VectorView& A){
  return MatrixViewMap(A.mdata+A.mrange.get_start(),
   A.nelem(), 1, StrideType(A.mrange.get_stride(), 1)); }

// Converts vector to eigen map row-view
MatrixViewMap MapToEigenRow(VectorView& A){return MapToEigen(A);}

// Converts vector to eigen map column-view
MatrixViewMap MapToEigenCol(VectorView& A){
  return MatrixViewMap(A.mdata+A.mrange.get_start(),
   1, A.nelem(), StrideType(1, A.mrange.get_stride())); }

// Special 4x4

// Converts matrix to eigen map
Matrix4x4ViewMap MapToEigen4x4(MatrixView& A){
  return Matrix4x4ViewMap(A.mdata+A.mrr.get_start()+A.mcr.get_start(),
                          4, 4, StrideType(A.mrr.get_stride(), A.mcr.get_stride())); }

// Converts constant matrix to constant eigen map
ConstMatrix4x4ViewMap MapToEigen4x4(const ConstMatrixView& A){
  return ConstMatrix4x4ViewMap(A.mdata+A.mrr.get_start()+A.mcr.get_start(),
                               4, 4, StrideType(A.mrr.get_stride(), A.mcr.get_stride())); }

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
Numeric debug_matrixview_get_elem (MatrixView& mv, Index r, Index c)
{
  return mv(r, c);
}

#endif
////////////////////////////////



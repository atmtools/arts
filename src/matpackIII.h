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
  Implementation of Tensors of Rank 3.

  The three dimensions are called: page, row, column.
  
  \author Stefan Buehler
  \date   2001-11-22
 */

#ifndef matpackIII_h
#define matpackIII_h

#include <iomanip>
#include "matpackI.h"

/** The outermost iterator class for rank 3 tensors. This takes into
    account the defined strided. */
class Iterator3D {
public:
  // Constructors:
  Iterator3D();
  Iterator3D(const Iterator3D& o);
  Iterator3D(const MatrixView& x, Index stride);

  // Operators:
  Iterator3D& operator++();
  bool operator!=(const Iterator3D& other) const;
  MatrixView* const operator->();
  MatrixView& operator*();
  
private:
  /** Current position. */
  MatrixView msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator3D. */
class ConstIterator3D {
public:
  // Constructors:
  ConstIterator3D();
  ConstIterator3D(const ConstIterator3D& o);
  ConstIterator3D(const ConstMatrixView& x, Index stride);

  // Operators:
  ConstIterator3D& operator++();
  bool operator!=(const ConstIterator3D& other) const;
  const ConstMatrixView* operator->() const;
  const ConstMatrixView& operator*() const;

private:
  /** Current position. */
  ConstMatrixView msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor3:
class Tensor3;


/** A constant view of a Tensor3.

    This, together with the derived class Tensor3View, contains the
    main implementation of a Tensor3. It defines the concepts of
    Tensor3View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor3View from a subrange of
    a Tensor3View.

    The three dimensions of the tensor are called: page, row, column.

    The class Tensor3 is just a special case of a Tensor3View
    which also allocates storage. */
class ConstTensor3View {
public:
  // Member functions:
  Index npages() const;
  Index nrows() const;
  Index ncols() const;

  // Const index operators:
  ConstTensor3View operator()( const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( Index p,        Index r,        const Range& c ) const;
  ConstVectorView  operator()( Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( const Range& p, Index r,        Index c        ) const;

  Numeric          operator()( Index p,        Index r,        Index c        ) const;

  // Functions returning iterators:
  ConstIterator3D begin() const;
  ConstIterator3D end() const;
  
  // Friends:
  friend class Tensor3View;


protected:
  // Constructors:
  ConstTensor3View();
  ConstTensor3View(Numeric *data,
                   const Range& p, const Range& r, const Range& c);
  ConstTensor3View(Numeric *data,
                   const Range& pp, const Range& pr, const Range& pc,
                   const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The page range of mdata that is actually used. */
  Range mpr;
  /** The row range of mdata that is actually used. */
  Range mrr;
  /** The column range of mdata that is actually used. */
  Range mcr;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The Tensor3View class

    This contains the main implementation of a Tensor3. It defines
    the concepts of Tensor3View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor3View
    from a subrange of a Tensor3View. 

    The class Tensor3 is just a special case of a Tensor3View
    which also allocates storage. */
class Tensor3View : public ConstTensor3View {
public:

  // Const index operators:
  ConstTensor3View operator()( const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( Index p,        Index r,        const Range& c ) const;
  ConstVectorView  operator()( Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( const Range& p, Index r,        Index c        ) const;

  Numeric          operator()( Index p,        Index r,        Index c        ) const;

  // Non-const index operators:

  Tensor3View operator()( const Range& p, const Range& r, const Range& c );

  MatrixView  operator()( const Range& p, const Range& r, Index c        );
  MatrixView  operator()( const Range& p, Index r,        const Range& c );
  MatrixView  operator()( Index p,        const Range& r, const Range& c );

  VectorView  operator()( Index p,        Index r,        const Range& c );
  VectorView  operator()( Index p,        const Range& r, Index c        );
  VectorView  operator()( const Range& p, Index r,        Index c        );

  Numeric&    operator()( Index p,        Index r,        Index c        );

  // Functions returning const iterators:
  ConstIterator3D begin() const;
  ConstIterator3D end() const;
  // Functions returning iterators:
  Iterator3D begin();
  Iterator3D end();
  
  // Assignment operators:
  Tensor3View& operator=(const ConstTensor3View& v);
  Tensor3View& operator=(const Tensor3View& v);
  Tensor3View& operator=(const Tensor3& v);
  Tensor3View& operator=(Numeric x);

  // Other operators:
  Tensor3View& operator*=(Numeric x);
  Tensor3View& operator/=(Numeric x);
  Tensor3View& operator+=(Numeric x);
  Tensor3View& operator-=(Numeric x);

  Tensor3View& operator*=(const ConstTensor3View& x);
  Tensor3View& operator/=(const ConstTensor3View& x);
  Tensor3View& operator+=(const ConstTensor3View& x);
  Tensor3View& operator-=(const ConstTensor3View& x);

  // Friends:
//   friend class VectorView;
//   friend ConstTensor3View transpose(ConstTensor3View m);
//   friend Tensor3View transpose(Tensor3View m);

protected:
  // Constructors:
  Tensor3View();
  Tensor3View(Numeric *data, const Range& p, const Range& r, const Range& c);
  Tensor3View(Numeric *data,
              const Range& pp, const Range& pr, const Range& pc,
              const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor3 class. This is a Tensor3View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor3View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor3 : public Tensor3View {
public:
  // Constructors:
  Tensor3();
  Tensor3(Index p, Index r, Index c);
  Tensor3(Index p, Index r, Index c, Numeric fill);
  Tensor3(const ConstTensor3View& v);
  Tensor3(const Tensor3& v);

  // Assignment operators:
  Tensor3& operator=(const Tensor3& x);
  Tensor3& operator=(Numeric x);

  // Resize function:
  void resize(Index p, Index r, Index c);

  // Destructor:
  ~Tensor3();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator3D origin,
                 const ConstIterator3D& end,
                 Iterator3D target);

inline void copy(Numeric x,
                 Iterator3D target,
                 const Iterator3D& end);



// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor3. */
typedef Array<Tensor3> ArrayOfTensor3;



// Functions for Iterator3D
// ------------------------

/** Default constructor. */
inline Iterator3D::Iterator3D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator3D::Iterator3D(const Iterator3D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator3D::Iterator3D(const MatrixView& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline Iterator3D& Iterator3D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool Iterator3D::operator!=(const Iterator3D& other) const
{
  if ( msv.mdata +
       msv.mrr.mstart +
       msv.mcr.mstart 
       !=
       other.msv.mdata +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline MatrixView* const Iterator3D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline MatrixView& Iterator3D::operator*()
{
  return msv;
}

// Functions for ConstIterator3D
// -----------------------------

/** Default constructor. */
inline ConstIterator3D::ConstIterator3D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator3D::ConstIterator3D(const ConstIterator3D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator3D::ConstIterator3D(const ConstMatrixView& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline ConstIterator3D& ConstIterator3D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool ConstIterator3D::operator!=(const ConstIterator3D& other) const
{
  if ( msv.mdata +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline const ConstMatrixView* ConstIterator3D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstMatrixView& ConstIterator3D::operator*() const
{
  return msv;
}



// Functions for ConstTensor3View:
// ------------------------------

/** Returns the number of pages. */
inline Index ConstTensor3View::npages() const
{
  return mpr.mextent;
}

/** Returns the number of rows. */
inline Index ConstTensor3View::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor3View::ncols() const
{
  return mcr.mextent;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior.  */
inline ConstTensor3View ConstTensor3View::operator()(const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  return ConstTensor3View(mdata, mpr, mrr, mcr, p, r, c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
inline ConstMatrixView ConstTensor3View::operator()(const Range& p,
                                                    const Range& r,
                                                    Index c) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata + mcr.mstart + c*mcr.mstride,
                         mpr, mrr, 
                         p,   r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
inline ConstMatrixView ConstTensor3View::operator()(const Range& p,
                                                    Index        r,
                                                    const Range& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return ConstMatrixView(mdata + mrr.mstart + r*mrr.mstride,
                         mpr, mcr, 
                         p,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
inline ConstMatrixView ConstTensor3View::operator()(Index p,
                                                    const Range& r,
                                                    const Range& c) const
{
  // Check that p is valid:
  assert( 0 <= p );
  assert( p < mpr.mextent );

  return ConstMatrixView(mdata + mpr.mstart + p*mpr.mstride,
                         mrr, mcr, 
                         r,   c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) */
inline ConstVectorView ConstTensor3View::operator()(Index p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that p and r are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstVectorView(mdata +
                         mpr.mstart + p*mpr.mstride +
                         mrr.mstart + r*mrr.mstride,
                         mcr, c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) */
inline ConstVectorView ConstTensor3View::operator()(Index p,
                                                    const Range& r,
                                                    Index c) const
{
  // Check that p and c are valid:
  assert( 0 <= p );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         mpr.mstart + p*mpr.mstride +
                         mcr.mstart + c*mcr.mstride,
                         mrr, r);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) */
inline ConstVectorView ConstTensor3View::operator()(const Range& p,
                                                    Index r,
                                                    Index c) const
{
  // Check that r and c are valid:
  assert( 0 <= r );
  assert( 0 <= c );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         mrr.mstart + r*mrr.mstride +
                         mcr.mstart + c*mcr.mstride,
                         mpr, p);
}

/** Plain const index operator. */
inline Numeric ConstTensor3View::operator()(Index p, Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=p );
  assert( 0<=r );
  assert( 0<=c );
  assert( p<mpr.mextent );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  return *( mdata +
            mpr.mstart + p*mpr.mstride +
            mrr.mstart + r*mrr.mstride +
            mcr.mstart + c*mcr.mstride );
}

/** Return const iterator to first page. */
inline ConstIterator3D ConstTensor3View::begin() const
{
  return ConstIterator3D(ConstMatrixView(mdata+mpr.mstart,
                                         mrr,
                                         mcr),
                         mpr.mstride);
}

/** Return const iterator behind last page. */
inline ConstIterator3D ConstTensor3View::end() const
{
  return ConstIterator3D( ConstMatrixView(mdata + mpr.mstart +
                                          (mpr.mextent)*mpr.mstride,
                                          mrr,
                                          mcr),
                          mpr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor3View::ConstTensor3View() :
  mpr(0,0,1), mrr(0,0,1), mcr(0,0,1), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor3 to initialize
    its own Tensor3View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
inline ConstTensor3View::ConstTensor3View(Numeric *data,
                                          const Range& pr,
                                          const Range& rr,
                                          const Range& cr) :
  mpr(pr),
  mrr(rr),
  mcr(cr),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
inline ConstTensor3View::ConstTensor3View(Numeric *data,
                                          const Range& pp,
                                          const Range& pr,
                                          const Range& pc,
                                          const Range& np,
                                          const Range& nr,
                                          const Range& nc) :
  mpr(pp,np),
  mrr(pr,nr),
  mcr(pc,nc),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Matrix to print each page in turn. */
inline std::ostream& operator<<(std::ostream& os, const ConstTensor3View& v)
{
  // Page iterators:
  ConstIterator3D ip=v.begin();
  const ConstIterator3D end_page=v.end();

  if ( ip!=end_page )
    {
      os << *ip;
      ++ip;
    }

  for ( ; ip!=end_page; ++ip )
    {
      os << "\n\n";
      os << *ip;
    }

  return os;
}


// Functions for Tensor3View:
// -------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline ConstTensor3View Tensor3View::operator()(const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor3View::operator()(const Range& p,
                                                    const Range& r,
                                                    Index c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor3View::operator()(const Range& p,
                                                    Index        r,
                                                    const Range& c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor3View::operator()(Index p,
                                                    const Range& r,
                                                    const Range& c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor3View::operator()(Index p,
                                                    Index r,
                                                    const Range& c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor3View::operator()(Index p,
                                                    const Range& r,
                                                    Index c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor3View::operator()(const Range& p,
                                                    Index r,
                                                    Index c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline Numeric Tensor3View::operator()(Index p, Index r, Index c) const
{
  return ConstTensor3View::operator()(p,r,c);  
}

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior.  */
inline Tensor3View Tensor3View::operator()(const Range& p,
                                           const Range& r,
                                           const Range& c) 
{
  return Tensor3View(mdata, mpr, mrr, mcr, p, r, c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
inline MatrixView Tensor3View::operator()(const Range& p,
                                          const Range& r,
                                          Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return MatrixView(mdata + mcr.mstart + c*mcr.mstride,
                    mpr, mrr, 
                    p,   r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
inline MatrixView Tensor3View::operator()(const Range& p,
                                          Index        r,
                                          const Range& c) 
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return MatrixView(mdata + mrr.mstart + r*mrr.mstride,
                    mpr, mcr, 
                    p,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
inline MatrixView Tensor3View::operator()(Index p,
                                          const Range& r,
                                          const Range& c) 
{
  // Check that p is valid:
  assert( 0 <= p );
  assert( p < mpr.mextent );

  return MatrixView(mdata + mpr.mstart + p*mpr.mstride,
                    mrr, mcr, 
                    r,   c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by two.) */
inline VectorView Tensor3View::operator()(Index p,
                                          Index r,
                                          const Range& c) 
{
  // Check that p and r are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return VectorView(mdata +
                    mpr.mstart + p*mpr.mstride +
                    mrr.mstart + r*mrr.mstride,
                    mcr, c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by two.) */
inline VectorView Tensor3View::operator()(Index p,
                                          const Range& r,
                                          Index c) 
{
  // Check that p and c are valid:
  assert( 0 <= p );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    mpr.mstart + p*mpr.mstride +
                    mcr.mstart + c*mcr.mstride,
                    mrr, r);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by two.) */
inline VectorView Tensor3View::operator()(const Range& p,
                                          Index r,
                                          Index c)
{
  // Check that r and r are valid:
  assert( 0 <= r );
  assert( 0 <= c );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    mrr.mstart + r*mrr.mstride +
                    mcr.mstart + c*mcr.mstride,
                    mpr, p);
}

/** Plain non-const index operator. */
inline Numeric& Tensor3View::operator()(Index p, Index r, Index c) 
{
  // Check if indices are valid:
  assert( 0<=p );
  assert( 0<=r );
  assert( 0<=c );
  assert( p<mpr.mextent );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  return *( mdata +
            mpr.mstart + p*mpr.mstride +
            mrr.mstart + r*mrr.mstride +
            mcr.mstart + c*mcr.mstride );
}

/** Return const iterator to first row. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator3D Tensor3View::begin() const
{
  return ConstTensor3View::begin();
}

/** Return const iterator behind last row. */
inline ConstIterator3D Tensor3View::end() const
{
  return ConstTensor3View::end();
}

/** Return iterator to first page. */
inline Iterator3D Tensor3View::begin()
{
  return Iterator3D(MatrixView(mdata+mpr.mstart,
                               mrr,
                               mcr),
                    mpr.mstride);
}

/** Return iterator behind last page. */
inline Iterator3D Tensor3View::end()
{
  return Iterator3D( MatrixView(mdata + mpr.mstart +
                                (mpr.mextent)*mpr.mstride,
                                mrr,
                                mcr),
                     mpr.mstride );
}

/** Assignment operator. This copies the data from another Tensor3View
    to this Tensor3View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor3View by
    setting its range. */
inline Tensor3View& Tensor3View::operator=(const ConstTensor3View& m)
{
  // Check that sizes are compatible:
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from Tensor3View to Tensor3View. This is a tricky
    one. The problem is that since Tensor3View is derived from
    ConstTensor3View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline Tensor3View& Tensor3View::operator=(const Tensor3View& m)
{
  // Check that sizes are compatible:
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Tensor3. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline Tensor3View& Tensor3View::operator=(const Tensor3& m)
{
  // Check that sizes are compatible:
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor3View will set all elements to this
    value. */
inline Tensor3View& Tensor3View::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

// Some little helper functions:
//------------------------------

inline Numeric add(Numeric x, Numeric y)
{
  return x+y;
}

/** Multiplication by scalar. */
inline Tensor3View& Tensor3View::operator*=(Numeric x)
{
  const Iterator3D ep=end();  
  for ( Iterator3D p=begin(); p!=ep ; ++p )
  {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
inline Tensor3View& Tensor3View::operator/=(Numeric x)
{
  const Iterator3D ep=end();  
  for ( Iterator3D p=begin(); p!=ep ; ++p )
  {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
inline Tensor3View& Tensor3View::operator+=(Numeric x)
{
  const Iterator3D ep=end();  
  for ( Iterator3D p=begin(); p!=ep ; ++p )
  {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
inline Tensor3View& Tensor3View::operator-=(Numeric x)
{
  const Iterator3D ep=end();  
  for ( Iterator3D p=begin(); p!=ep ; ++p )
  {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor3. */
inline Tensor3View& Tensor3View::operator*=(const ConstTensor3View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator3D  xp = x.begin();
  Iterator3D        p = begin();
  const Iterator3D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p *= *xp;
    }
  return *this;
}

/** Element-vise division by another Tensor3. */
inline Tensor3View& Tensor3View::operator/=(const ConstTensor3View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator3D  xp = x.begin();
  Iterator3D        p = begin();
  const Iterator3D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p /= *xp;
    }
  return *this;
}

/** Element-vise addition of another Tensor3. */
inline Tensor3View& Tensor3View::operator+=(const ConstTensor3View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator3D  xp = x.begin();
  Iterator3D        p = begin();
  const Iterator3D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p += *xp;
    }
  return *this;
}

/** Element-vise subtraction of another Tensor3. */
inline Tensor3View& Tensor3View::operator-=(const ConstTensor3View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator3D  xp = x.begin();
  Iterator3D        p = begin();
  const Iterator3D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p -= *xp;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Tensor3. */
inline Tensor3View::Tensor3View() :
  ConstTensor3View()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor3 to initialize its
    own Tensor3View part. The row range rr must have a
    stride to account for the length of one row. */
inline Tensor3View::Tensor3View(Numeric *data,
                                const Range& pr,
                                const Range& rr,
                                const Range& cr) :
  ConstTensor3View(data, pr, rr, cr)
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
    \param p Previous range.
    \param n New Range.  */
inline Tensor3View::Tensor3View(Numeric *data,
                                const Range& pp,
                                const Range& pr,
                                const Range& pc,
                                const Range& np,
                                const Range& nr,
                                const Range& nc) :
  ConstTensor3View(data,pp,pr,pc,np,nr,nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
inline void copy(ConstIterator3D origin,
                 const ConstIterator3D& end,
                 Iterator3D target)
{
  for ( ; origin!=end ; ++origin,++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy(origin->begin(),
           origin->end(),
           target->begin());
    }
}

/** Copy a scalar to all elements. */
inline void copy(Numeric x,
                 Iterator3D target,
                 const Iterator3D& end)
{
  for ( ; target!=end ; ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy(x,target->begin(),target->end());
    }
}


// Functions for Tensor3:
// ---------------------

/** Default constructor. */
inline Tensor3::Tensor3() :
  Tensor3View::Tensor3View()
{
  // Nothing to do here. However, note that the default constructor
  // for Tensor3View has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized. 
}

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
inline Tensor3::Tensor3(Index p, Index r, Index c) :
  Tensor3View( new Numeric[p*r*c],
             Range(0,p,r*c),
             Range(0,r,c),
             Range(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Tensor3::Tensor3(Index p, Index r, Index c, Numeric fill) :
  Tensor3View( new Numeric[p*r*c],
              Range(0,p,r*c),
              Range(0,r,c),
              Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x=mdata; x<mdata+p*r*c; ++x )
    *x = fill;
}

/** Copy constructor from Tensor3View. This automatically sets the size
    and copies the data. */
inline Tensor3::Tensor3(const ConstTensor3View& m) :
  Tensor3View( new Numeric[m.npages()*m.nrows()*m.ncols()],
             Range( 0, m.npages(), m.nrows()*m.ncols() ),
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Tensor3. This automatically sets the size
    and copies the data. */
inline Tensor3::Tensor3(const Tensor3& m) :
  Tensor3View( new Numeric[m.npages()*m.nrows()*m.ncols()],
             Range( 0, m.npages(), m.nrows()*m.ncols() ),
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty tensor, then it will have
  // dimensions of size 0. But these are used to initialize the stride
  // for higher dimensions! Thus, this method has to be consistent
  // with the behaviour of Range::Range. For now, Range::Range allows
  // also stride 0.
  copy(m.begin(),m.end(),begin());
}

/** Assignment operator from another tensor. It is important that this
    operator exists. Otherwise the = operator seems to copy references
    instead of content in some cases. 

    The Behavior of this one is a bit special: If the size of the
    target tensor is 0 then it will be automatically resized to match
    (this is needed to have the correct initialization for constructed
    classes that use the assignment operator to initialize their
    data). 
*/
inline Tensor3& Tensor3::operator=(const Tensor3& m)
{
  //  cout << "Tensor3 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "             n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid tensor, so we just
  // have to check one.
  if ( 0 == mcr.mextent )
    {
      // Adjust if previously empty.
      resize( m.mpr.mextent, m.mrr.mextent, m.mcr.mextent ); 
    }
  else
    {
      // Check that sizes are compatible:
      assert( mpr.mextent==m.mpr.mextent );
      assert( mrr.mextent==m.mrr.mextent );
      assert( mcr.mextent==m.mcr.mextent );
    }

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
inline Tensor3& Tensor3::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
inline void Tensor3::resize(Index p, Index r, Index c)
{
  assert( 0<=p );
  assert( 0<=r );
  assert( 0<=c );

  if ( mpr.mextent!=p ||
       mrr.mextent!=r ||
       mcr.mextent!=c )
    {
      delete mdata;
      mdata = new Numeric[p*r*c];

      mpr.mstart = 0;
      mpr.mextent = p;
      mpr.mstride = r*c;

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
}

/** Destructor for Tensor3. This is important, since Tensor3 uses new to
    allocate storage. */
inline Tensor3::~Tensor3()
{
//   cout << "Destroying a Tensor3:\n"
//        << *this << "\n........................................\n";
  delete mdata;
}


/** A generic transform function for tensors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for tensors! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    The two views may be the same one, in which case the
    conversion happens in place. 

    \retval   y   the results of the function acting on each element of x
    \param    my_func a function (e.g., sqrt)
    \param    x   a tensor */
inline void transform( Tensor3View y,
                       double (&my_func)(double),
                       ConstTensor3View x )
{
  // Check dimensions:
  assert( y.npages() == x.npages() );
  assert( y.nrows()  == x.nrows()  );
  assert( y.ncols()  == x.ncols()  );

  const ConstIterator3D xe = x.end();
  ConstIterator3D       xi = x.begin();
  Iterator3D            yi = y.begin();
  for ( ; xi!=xe; ++xi, ++yi )
    {
      // Use the transform function of lower dimensional tensors
      // recursively:
      transform(*yi,my_func,*xi);
    }
}

/** Max function, tensor version. */
inline Numeric max(const ConstTensor3View& x)
{
  const ConstIterator3D xe = x.end();
  ConstIterator3D       xi = x.begin();

  // Initial value for max:
  Numeric themax = max(*xi);
  ++xi;

  for ( ; xi!=xe ; ++xi )
    {
      // Use the max function of lower dimensional tensors
      // recursively:
      Numeric maxi = max(*xi);
      if ( maxi > themax )
        themax = maxi;
    }

  return themax;
}

/** Min function, tensor version. */
inline Numeric min(const ConstTensor3View& x)
{
  const ConstIterator3D xe = x.end();
  ConstIterator3D       xi = x.begin();

  // Initial value for min:
  Numeric themin = min(*xi);
  ++xi;

  for ( ; xi!=xe ; ++xi )
    {
      // Use the min function of lower dimensional tensors
      // recursively:
      Numeric mini = min(*xi);
      if ( mini < themin )
        themin = mini;
    }

  return themin;
}



#endif    // matpackIII_h

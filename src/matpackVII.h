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
  Implementation of Tensors of Rank 7.

  Dimensions are called: library, vitrine, shelf, book, page, row, column.
  or short:              l,       v,       s,     b,    p,    r,   c
  
  \author Stefan Buehler
  \date   2001-11-22
 */

#ifndef matpackVII_h
#define matpackVII_h

#include <iomanip>
#include "matpackI.h"

/** The outermost iterator class for rank 7 tensors. This takes into
    account the defined strided. */
class Iterator7D {
public:
  // Constructors:
  Iterator7D();
  Iterator7D(const Iterator7D& o);
  Iterator7D(const Tensor6View& x, Index stride);

  // Operators:
  Iterator7D& operator++();
  bool operator!=(const Iterator7D& other) const;
  MatrixView* const operator->();
  MatrixView& operator*();
  
private:
  /** Current position. */
  Tensor6View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator7D. */
class ConstIterator7D {
public:
  // Constructors:
  ConstIterator7D();
  ConstIterator7D(const ConstIterator7D& o);
  ConstIterator7D(const ConstTensor6View& x, Index stride);

  // Operators:
  ConstIterator7D& operator++();
  bool operator!=(const ConstIterator7D& other) const;
  const ConstMatrixView* operator->() const;
  const ConstMatrixView& operator*() const;

private:
  /** Current position. */
  ConstMatrixView msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor7:
class Tensor7;


/** A constant view of a Tensor7.

    This, together with the derived class Tensor7View, contains the
    main implementation of a Tensor7. It defines the concepts of
    Tensor7View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor7View from a subrange of
    a Tensor7View.

    Dimensions are called: library, vitrine, shelf, book, page, row, column.
    or short:              l,       v,       s,     b,    p,    r,   c

    The class Tensor7 is just a special case of a Tensor7View
    which also allocates storage. */
class ConstTensor7View {
public:
  // Member functions:
  Index nlibraries() const;
  Index nvitrines()  const;
  Index nshelves()   const;
  Index nbooks()     const;
  Index npages()     const;
  Index nrows()      const;
  Index ncols()      const;

  // Const index operators:
  ConstTensor7View operator()( const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( Index p,        Index r,        const Range& c ) const;
  ConstVectorView  operator()( Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( const Range& p, Index r,        Index c        ) const;

  Numeric          operator()( Index p,        Index r,        Index c        ) const;

  // Functions returning iterators:
  ConstIterator7D begin() const;
  ConstIterator7D end() const;
  
  // Friends:
  friend class Tensor7View;
  friend class ConstIterator4D;
  friend class ConstTensor4View;
  friend class ConstTensor5View;


protected:
  // Constructors:
  ConstTensor7View();
  ConstTensor7View(Numeric *data,
                   const Range& p, const Range& r, const Range& c);
  ConstTensor7View(Numeric *data,
                   const Range& pp, const Range& pr, const Range& pc,
                   const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The library range of mdata that is actually used. */
  Range mlr;
  /** The vitrine range of mdata that is actually used. */
  Range mvr;
  /** The shelf range of mdata that is actually used. */
  Range msr;
  /** The book range of mdata that is actually used. */
  Range mbr;
  /** The page range of mdata that is actually used. */
  Range mpr;
  /** The row range of mdata that is actually used. */
  Range mrr;
  /** The column range of mdata that is actually used. */
  Range mcr;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The Tensor7View class

    This contains the main implementation of a Tensor7. It defines
    the concepts of Tensor7View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor7View
    from a subrange of a Tensor7View. 

    The class Tensor7 is just a special case of a Tensor7View
    which also allocates storage. */
class Tensor7View : public ConstTensor7View {
public:

  // Const index operators:
  ConstTensor7View operator()( const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( Index p,        Index r,        const Range& c ) const;
  ConstVectorView  operator()( Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( const Range& p, Index r,        Index c        ) const;

  Numeric          operator()( Index p,        Index r,        Index c        ) const;

  // Non-const index operators:

  Tensor7View operator()( const Range& p, const Range& r, const Range& c );

  MatrixView  operator()( const Range& p, const Range& r, Index c        );
  MatrixView  operator()( const Range& p, Index r,        const Range& c );
  MatrixView  operator()( Index p,        const Range& r, const Range& c );

  VectorView  operator()( Index p,        Index r,        const Range& c );
  VectorView  operator()( Index p,        const Range& r, Index c        );
  VectorView  operator()( const Range& p, Index r,        Index c        );

  Numeric&    operator()( Index p,        Index r,        Index c        );

  // Functions returning const iterators:
  ConstIterator7D begin() const;
  ConstIterator7D end() const;
  // Functions returning iterators:
  Iterator7D begin();
  Iterator7D end();
  
  // Assignment operators:
  Tensor7View& operator=(const ConstTensor7View& v);
  Tensor7View& operator=(const Tensor7View& v);
  Tensor7View& operator=(const Tensor7& v);
  Tensor7View& operator=(Numeric x);

  // Other operators:
  Tensor7View& operator*=(Numeric x);
  Tensor7View& operator/=(Numeric x);
  Tensor7View& operator+=(Numeric x);
  Tensor7View& operator-=(Numeric x);

  Tensor7View& operator*=(const ConstTensor7View& x);
  Tensor7View& operator/=(const ConstTensor7View& x);
  Tensor7View& operator+=(const ConstTensor7View& x);
  Tensor7View& operator-=(const ConstTensor7View& x);

  // Friends:
//   friend class VectorView;
//   friend ConstTensor7View transpose(ConstTensor7View m);
//   friend Tensor7View transpose(Tensor7View m);
  friend class Iterator4D;
  friend class Tensor4View;
  friend class Tensor5View;

protected:
  // Constructors:
  Tensor7View();
  Tensor7View(Numeric *data, const Range& p, const Range& r, const Range& c);
  Tensor7View(Numeric *data,
	      const Range& pp, const Range& pr, const Range& pc,
	      const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor7 class. This is a Tensor7View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor7View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor7 : public Tensor7View {
public:
  // Constructors:
  Tensor7();
  Tensor7(Index p, Index r, Index c);
  Tensor7(Index p, Index r, Index c, Numeric fill);
  Tensor7(const ConstTensor7View& v);
  Tensor7(const Tensor7& v);

  // Assignment operators:
  Tensor7& operator=(const Tensor7& x);
  Tensor7& operator=(Numeric x);

  // Resize function:
  void resize(Index p, Index r, Index c);

  // Destructor:
  ~Tensor7();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator7D origin,
                 const ConstIterator7D& end,
                 Iterator7D target);

inline void copy(Numeric x,
                 Iterator7D target,
                 const Iterator7D& end);



// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor7. */
typedef Array<Tensor7> ArrayOfTensor7;



// Functions for Iterator7D
// ------------------------

/** Default constructor. */
inline Iterator7D::Iterator7D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator7D::Iterator7D(const Iterator7D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator7D::Iterator7D(const Tensor6View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline Iterator7D& Iterator7D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. 
    FIXME: Is it really necessary to have such a complicated check
    here? It could be sufficient to just test
    msv.mdata!=other.msv.mdata. */
inline bool Iterator7D::operator!=(const Iterator7D& other) const
{
  if ( msv.mdata +
       msv.mvr.mstart +
       msv.msr.mstart +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart 
       !=
       other.msv.mdata +
       other.msv.mvr.mstart +
       other.msv.msr.mstart +
       other.msv.mbr.mstart +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline MatrixView* const Iterator7D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline MatrixView& Iterator7D::operator*()
{
  return msv;
}

// Functions for ConstIterator7D
// -----------------------------

/** Default constructor. */
inline ConstIterator7D::ConstIterator7D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator7D::ConstIterator7D(const ConstIterator7D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator7D::ConstIterator7D(const ConstTensor6View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline ConstIterator7D& ConstIterator7D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. 
    FIXME: Is it really necessary to have such a complicated check
    here? It could be sufficient to just test
    msv.mdata!=other.msv.mdata. */
inline bool ConstIterator7D::operator!=(const ConstIterator7D& other) const
{
  if ( msv.mdata +
       msv.mvr.mstart +
       msv.msr.mstart +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mvr.mstart +
       other.msv.msr.mstart +
       other.msv.mbr.mstart +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline const ConstMatrixView* ConstIterator7D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstMatrixView& ConstIterator7D::operator*() const
{
  return msv;
}



// Functions for ConstTensor7View:
// ------------------------------

/** Returns the number of libraries. */
inline Index ConstTensor7View::nlibraries() const
{
  return mlr.mextent;
}

/** Returns the number of vitrines. */
inline Index ConstTensor7View::nvitrines() const
{
  return mvr.mextent;
}

/** Returns the number of shelves. */
inline Index ConstTensor7View::nshelves() const
{
  return msr.mextent;
}

/** Returns the number of books. */
inline Index ConstTensor7View::nbooks() const
{
  return mbr.mextent;
}

/** Returns the number of pages. */
inline Index ConstTensor7View::npages() const
{
  return mpr.mextent;
}

/** Returns the number of rows. */
inline Index ConstTensor7View::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor7View::ncols() const
{
  return mcr.mextent;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor7. This allows
    correct recursive behavior.  */
inline ConstTensor7View ConstTensor7View::operator()(const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  return ConstTensor7View(mdata, mpr, mrr, mcr, p, r, c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
inline ConstMatrixView ConstTensor7View::operator()(const Range& p,
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
inline ConstMatrixView ConstTensor7View::operator()(const Range& p,
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
inline ConstMatrixView ConstTensor7View::operator()(Index p,
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
inline ConstVectorView ConstTensor7View::operator()(Index p,
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
inline ConstVectorView ConstTensor7View::operator()(Index p,
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
inline ConstVectorView ConstTensor7View::operator()(const Range& p,
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
inline Numeric ConstTensor7View::operator()(Index p, Index r, Index c) const
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
inline ConstIterator7D ConstTensor7View::begin() const
{
  return ConstIterator7D(ConstMatrixView(mdata+mpr.mstart,
                                         mrr,
                                         mcr),
                         mpr.mstride);
}

/** Return const iterator behind last page. */
inline ConstIterator7D ConstTensor7View::end() const
{
  return ConstIterator7D( ConstMatrixView(mdata + mpr.mstart +
                                          (mpr.mextent)*mpr.mstride,
                                          mrr,
                                          mcr),
                          mpr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor7View::ConstTensor7View() :
  mpr(0,0,1), mrr(0,0,1), mcr(0,0,1), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor7 to initialize
    its own Tensor7View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
inline ConstTensor7View::ConstTensor7View(Numeric *data,
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
inline ConstTensor7View::ConstTensor7View(Numeric *data,
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
inline std::ostream& operator<<(std::ostream& os, const ConstTensor7View& v)
{
  // Page iterators:
  ConstIterator7D ip=v.begin();
  const ConstIterator7D end_page=v.end();

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


// Functions for Tensor7View:
// -------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor7. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline ConstTensor7View Tensor7View::operator()(const Range& p,
						const Range& r,
						const Range& c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor7View::operator()(const Range& p,
                                                    const Range& r,
                                                    Index c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor7View::operator()(const Range& p,
                                                    Index        r,
                                                    const Range& c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor7View::operator()(Index p,
                                                    const Range& r,
                                                    const Range& c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor7View::operator()(Index p,
                                                    Index r,
                                                    const Range& c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor7View::operator()(Index p,
                                                    const Range& r,
                                                    Index c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor7View::operator()(const Range& p,
                                                    Index r,
                                                    Index c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline Numeric Tensor7View::operator()(Index p, Index r, Index c) const
{
  return ConstTensor7View::operator()(p,r,c);  
}

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor7. This allows
    correct recursive behavior.  */
inline Tensor7View Tensor7View::operator()(const Range& p,
					   const Range& r,
					   const Range& c) 
{
  return Tensor7View(mdata, mpr, mrr, mcr, p, r, c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
inline MatrixView Tensor7View::operator()(const Range& p,
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
inline MatrixView Tensor7View::operator()(const Range& p,
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
inline MatrixView Tensor7View::operator()(Index p,
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
inline VectorView Tensor7View::operator()(Index p,
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
inline VectorView Tensor7View::operator()(Index p,
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
inline VectorView Tensor7View::operator()(const Range& p,
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
inline Numeric& Tensor7View::operator()(Index p, Index r, Index c) 
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
inline ConstIterator7D Tensor7View::begin() const
{
  return ConstTensor7View::begin();
}

/** Return const iterator behind last row. */
inline ConstIterator7D Tensor7View::end() const
{
  return ConstTensor7View::end();
}

/** Return iterator to first page. */
inline Iterator7D Tensor7View::begin()
{
  return Iterator7D(MatrixView(mdata+mpr.mstart,
			       mrr,
			       mcr),
		    mpr.mstride);
}

/** Return iterator behind last page. */
inline Iterator7D Tensor7View::end()
{
  return Iterator7D( MatrixView(mdata + mpr.mstart +
				(mpr.mextent)*mpr.mstride,
				mrr,
				mcr),
		     mpr.mstride );
}

/** Assignment operator. This copies the data from another Tensor7View
    to this Tensor7View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor7View by
    setting its range. */
inline Tensor7View& Tensor7View::operator=(const ConstTensor7View& m)
{
  // Check that sizes are compatible:
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from Tensor7View to Tensor7View. This is a tricky
    one. The problem is that since Tensor7View is derived from
    ConstTensor7View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline Tensor7View& Tensor7View::operator=(const Tensor7View& m)
{
  // Check that sizes are compatible:
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Tensor7. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline Tensor7View& Tensor7View::operator=(const Tensor7& m)
{
  // Check that sizes are compatible:
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor7View will set all elements to this
    value. */
inline Tensor7View& Tensor7View::operator=(Numeric x)
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
inline Tensor7View& Tensor7View::operator*=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
inline Tensor7View& Tensor7View::operator/=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
inline Tensor7View& Tensor7View::operator+=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
inline Tensor7View& Tensor7View::operator-=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor7. */
inline Tensor7View& Tensor7View::operator*=(const ConstTensor7View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()	 );
  assert( ncols()  == x.ncols()	 );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p *= *xp;
    }
  return *this;
}

/** Element-vise division by another Tensor7. */
inline Tensor7View& Tensor7View::operator/=(const ConstTensor7View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()	 );
  assert( ncols()  == x.ncols()	 );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p /= *xp;
    }
  return *this;
}

/** Element-vise addition of another Tensor7. */
inline Tensor7View& Tensor7View::operator+=(const ConstTensor7View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()	 );
  assert( ncols()  == x.ncols()	 );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p += *xp;
    }
  return *this;
}

/** Element-vise subtraction of another Tensor7. */
inline Tensor7View& Tensor7View::operator-=(const ConstTensor7View& x)
{
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()	 );
  assert( ncols()  == x.ncols()	 );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p -= *xp;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Tensor7. */
inline Tensor7View::Tensor7View() :
  ConstTensor7View()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor7 to initialize its
    own Tensor7View part. The row range rr must have a
    stride to account for the length of one row. */
inline Tensor7View::Tensor7View(Numeric *data,
				const Range& pr,
				const Range& rr,
				const Range& cr) :
  ConstTensor7View(data, pr, rr, cr)
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
inline Tensor7View::Tensor7View(Numeric *data,
				const Range& pp,
				const Range& pr,
				const Range& pc,
				const Range& np,
				const Range& nr,
				const Range& nc) :
  ConstTensor7View(data,pp,pr,pc,np,nr,nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
inline void copy(ConstIterator7D origin,
                 const ConstIterator7D& end,
                 Iterator7D target)
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
                 Iterator7D target,
                 const Iterator7D& end)
{
  for ( ; target!=end ; ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy(x,target->begin(),target->end());
    }
}


// Functions for Tensor7:
// ---------------------

/** Default constructor. */
inline Tensor7::Tensor7() :
  Tensor7View::Tensor7View()
{
  // Nothing to do here. However, note that the default constructor
  // for Tensor7View has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized. 
}

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
inline Tensor7::Tensor7(Index p, Index r, Index c) :
  Tensor7View( new Numeric[p*r*c],
             Range(0,p,r*c),
             Range(0,r,c),
             Range(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Tensor7::Tensor7(Index p, Index r, Index c, Numeric fill) :
  Tensor7View( new Numeric[p*r*c],
              Range(0,p,r*c),
              Range(0,r,c),
              Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x=mdata; x<mdata+p*r*c; ++x )
    *x = fill;
}

/** Copy constructor from Tensor7View. This automatically sets the size
    and copies the data. */
inline Tensor7::Tensor7(const ConstTensor7View& m) :
  Tensor7View( new Numeric[m.npages()*m.nrows()*m.ncols()],
             Range( 0, m.npages(), m.nrows()*m.ncols() ),
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Tensor7. This automatically sets the size
    and copies the data. */
inline Tensor7::Tensor7(const Tensor7& m) :
  Tensor7View( new Numeric[m.npages()*m.nrows()*m.ncols()],
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
inline Tensor7& Tensor7::operator=(const Tensor7& m)
{
  //  cout << "Tensor7 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
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
inline Tensor7& Tensor7::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
inline void Tensor7::resize(Index p, Index r, Index c)
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

/** Destructor for Tensor7. This is important, since Tensor7 uses new to
    allocate storage. */
inline Tensor7::~Tensor7()
{
//   cout << "Destroying a Tensor7:\n"
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
inline void transform( Tensor7View y,
                       double (&my_func)(double),
                       ConstTensor7View x )
{
  // Check dimensions:
  assert( y.npages() ==	x.npages() );
  assert( y.nrows()  ==	x.nrows()  );
  assert( y.ncols()  ==	x.ncols()  );

  const ConstIterator7D xe = x.end();
  ConstIterator7D       xi = x.begin();
  Iterator7D            yi = y.begin();
  for ( ; xi!=xe; ++xi, ++yi )
    {
      // Use the transform function of lower dimensional tensors
      // recursively:
      transform(*yi,my_func,*xi);
    }
}

/** Max function, tensor version. */
inline Numeric max(const ConstTensor7View& x)
{
  const ConstIterator7D xe = x.end();
  ConstIterator7D       xi = x.begin();

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
inline Numeric min(const ConstTensor7View& x)
{
  const ConstIterator7D xe = x.end();
  ConstIterator7D       xi = x.begin();

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



#endif    // matpackVII_h

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
  Implementation of Tensors of Rank 3 to 5.
  
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

    The class Tensor3 is just a special case of a Tensor3View
    which also allocates storage. */
class ConstTensor3View {
public:
  // Member functions:
  Index nelem1() const;
  Index nelem2() const;
  Index nelem3() const;

  // Const index operators:
  Numeric  operator()(Index r, Index c) const;
  ConstTensor3View operator()(const Range& r, const Range& c) const;
  ConstVectorView operator()(const Range& r, Index c) const;
  ConstVectorView operator()(Index r, const Range& c) const;

  // Functions returning iterators:
  ConstIterator2D begin() const;
  ConstIterator2D end() const;
  
  // Friends:
  friend class ConstVectorView;
  friend class Tensor3View;
  friend ConstTensor3View transpose(ConstTensor3View m);

protected:
  // Constructors:
  ConstTensor3View();
  ConstTensor3View(Numeric *data, const Range& r, const Range& c);
  ConstTensor3View(Numeric *data,
		  const Range& pr, const Range& pc,
		  const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The page range of mdata that is actually used. (first dimension)*/
  //  FIXME: I was here. Should I use names like page, row, column, or just number the dimensions? 
  Range mrr;
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
  Numeric  operator()(Index r, Index c) const;
  ConstTensor3View operator()(const Range& r, const Range& c) const;
  ConstVectorView operator()(const Range& r, Index c) const;
  ConstVectorView operator()(Index r, const Range& c) const;
  // Index Operators:
  Numeric& operator()(Index r, Index c);
  Tensor3View operator()(const Range& r, const Range& c);
  VectorView operator()(const Range& r, Index c);
  VectorView operator()(Index r, const Range& c);

  // Functions returning const iterators:
  ConstIterator2D begin() const;
  ConstIterator2D end() const;
  // Functions returning iterators:
  Iterator2D begin();
  Iterator2D end();
  
  // Assignment operators:
  Tensor3View& operator=(const ConstTensor3View& v);
  Tensor3View& operator=(const Tensor3View& v);
  Tensor3View& operator=(const Tensor3& v);
  Tensor3View& operator=(const ConstVectorView& v);
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

  Tensor3View& operator*=(const ConstVectorView& x);
  Tensor3View& operator/=(const ConstVectorView& x);
  Tensor3View& operator+=(const ConstVectorView& x);
  Tensor3View& operator-=(const ConstVectorView& x);

  // Friends:
  friend class VectorView;
  friend ConstTensor3View transpose(ConstTensor3View m);
  friend Tensor3View transpose(Tensor3View m);

protected:
  // Constructors:
  Tensor3View();
  Tensor3View(Numeric *data, const Range& r, const Range& c);
  Tensor3View(Numeric *data,
	     const Range& pr, const Range& pc,
	     const Range& nr, const Range& nc);
};

/** The Tensor3 class. This is a Tensor3View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor3View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operator from scalar.
    3. Resize function. */
class Tensor3 : public Tensor3View {
public:
  // Constructors:
  Tensor3();
  Tensor3(Index r, Index c);
  Tensor3(Index r, Index c, Numeric fill);
  Tensor3(const ConstTensor3View& v);
  Tensor3(const Tensor3& v);

  // Assignment operators:
  Tensor3& operator=(const Tensor3& x);
  Tensor3& operator=(Numeric x);
  Tensor3& operator=(const ConstVectorView& v);

  // Resize function:
  void resize(Index r, Index c);

  // Destructor:
  ~Tensor3();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator1D origin,
		 const ConstIterator1D& end,
		 Iterator1D target);

inline void copy(Numeric x,
		 Iterator1D target,
		 const Iterator1D& end);

inline void copy(ConstIterator2D origin,
		 const ConstIterator2D& end,
		 Iterator2D target);

inline void copy(Numeric x,
		 Iterator2D target,
		 const Iterator2D& end);



// Declare the existance of class Array:
template<class base>
class Array;

/** An array of vectors. */
typedef Array<Vector> ArrayOfVector;

/** An array of matrices. */
typedef Array<Matrix> ArrayOfMatrix;


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
  if ( msv.mdata + msv.mrange.mstart !=
       other.msv.mdata + other.msv.mrange.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline VectorView* const Iterator3D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline VectorView& Iterator3D::operator*()
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
  if ( msv.mdata + msv.mrange.mstart !=
       other.msv.mdata + other.msv.mrange.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline const ConstVectorView* ConstIterator3D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstVectorView& ConstIterator3D::operator*() const
{
  return msv;
}



// Functions for ConstTensor3View:
// ------------------------------

// FIXME: I was here. Should I use names like page, row, column, or just number the dimensions? 

/** Returns the number of rows. */
inline Index ConstTensor3View::nelem1() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor3View::ncols() const
{
  return mcr.mextent;
}

/** Plain const index operator. */
inline Numeric ConstTensor3View::operator()(Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  return *( mdata +
	    mrr.mstart +
	    r*mrr.mstride +
	    mcr.mstart +
	    c*mcr.mstride );
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior.  */
inline ConstTensor3View ConstTensor3View::operator()(const Range& r,
						   const Range& c) const
{
  return ConstTensor3View(mdata, mrr, mcr, r, c);
}

/** Const index operator returning a column as an object of type
    ConstVectorView.

    \param r A range of rows.
    \param c Index of selected column */
inline ConstVectorView ConstTensor3View::operator()(const Range& r, Index c) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c <= mcr.mstart+(mcr.mextent-1)*mcr.mstride );

  return ConstVectorView(mdata + mcr.mstart + c*mcr.mstride,
			 mrr, r);
}

/** Const index operator returning a row as an object of type
    ConstVectorView.

    \param r Index of selected row.
    \param c Range of columns */
inline ConstVectorView ConstTensor3View::operator()(Index r, const Range& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r <  mrr.mstart+(mrr.mextent-1)*mrr.mstride );

  return ConstVectorView(mdata + mrr.mstart + r*mrr.mstride,
			 mcr, c);
}

/** Return const iterator to first row. */
inline ConstIterator2D ConstTensor3View::begin() const
{
  return ConstIterator2D(ConstVectorView(mdata+mrr.mstart,
					 mcr),
			 mrr.mstride);
}

/** Return const iterator behind last row. */
inline ConstIterator2D ConstTensor3View::end() const
{
  return ConstIterator2D( ConstVectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
					  mcr),
			  mrr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor3View::ConstTensor3View() :
  mrr(0,0,1), mcr(0,0,1), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor3 to initialize its
    own Tensor3View part. The row range rr must have a
    stride to account for the length of one row. */
inline ConstTensor3View::ConstTensor3View(Numeric *data,
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
    \param p Previous range.
    \param n New Range.  */
inline ConstTensor3View::ConstTensor3View(Numeric *data,
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
inline std::ostream& operator<<(std::ostream& os, const ConstTensor3View& v)
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


// Functions for Tensor3View:
// -------------------------

/** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline Numeric Tensor3View::operator()(Index r, Index c) const
{
  return ConstTensor3View::operator()(r,c);
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline ConstTensor3View Tensor3View::operator()(const Range& r, const Range& c) const
{
  return ConstTensor3View::operator()(r,c);  
}

/** Const index operator returning a column as an object of type
    ConstVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r A range of rows.
    \param c Index of selected column */
inline ConstVectorView Tensor3View::operator()(const Range& r, Index c) const
{
  return ConstTensor3View::operator()(r,c);
}

/** Const index operator returning a row as an object of type
    ConstVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r Index of selected row.
    \param c Range of columns */
inline ConstVectorView Tensor3View::operator()(Index r, const Range& c) const
{
  return ConstTensor3View::operator()(r,c);
}

/** Plain index operator. */
inline Numeric& Tensor3View::operator()(Index r, Index c)
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  return *( mdata +
	    mrr.mstart +
	    r*mrr.mstride +
	    mcr.mstart +
	    c*mcr.mstride );
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Tensor3. This allows correct
    recursive behavior.  */
inline Tensor3View Tensor3View::operator()(const Range& r, const Range& c)
{
  return Tensor3View(mdata, mrr, mcr, r, c);
}

/** Index operator returning a column as an object of type VectorView.

    \param r A range of rows.
    \param c Index of selected column */
inline VectorView Tensor3View::operator()(const Range& r, Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c <= mcr.mstart+(mcr.mextent-1)*mcr.mstride );

  return VectorView(mdata + mcr.mstart + c*mcr.mstride,
		    mrr, r);
}

/** Index operator returning a row as an object of type VectorView.

    \param r Index of selected row.
    \param c Range of columns */
inline VectorView Tensor3View::operator()(Index r, const Range& c)
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r <  mrr.mstart+(mrr.mextent-1)*mrr.mstride );

  return VectorView(mdata + mrr.mstart + r*mrr.mstride,
		    mcr, c);
}

/** Return const iterator to first row. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator2D Tensor3View::begin() const
{
  return ConstTensor3View::begin();
}

/** Return const iterator behind last row. */
inline ConstIterator2D Tensor3View::end() const
{
  return ConstTensor3View::end();
}

/** Return iterator to first row. */
inline Iterator2D Tensor3View::begin()
{
  return Iterator2D(VectorView(mdata+mrr.mstart, mcr),
		    mrr.mstride);
}

/** Return iterator behind last row. */
inline Iterator2D Tensor3View::end()
{
  return Iterator2D( VectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
				mcr),
		     mrr.mstride );
}

/** Assignment operator. This copies the data from another Tensor3View
    to this Tensor3View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor3View by
    setting its range. */
inline Tensor3View& Tensor3View::operator=(const ConstTensor3View& m)
{
  // Check that sizes are compatible:
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
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a vector. This copies the data from a VectorView
    to this Tensor3View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor3View by
    setting its range. */
inline Tensor3View& Tensor3View::operator=(const ConstVectorView& v)
{
  // Check that sizes are compatible:
  assert( mrr.mextent==v.nelem() );
  assert( mcr.mextent==1         );
  //  dummy = ConstTensor3View(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstTensor3View dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor3View will set all elements to this
    value. */
inline Tensor3View& Tensor3View::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
inline Tensor3View& Tensor3View::operator*=(Numeric x)
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
inline Tensor3View& Tensor3View::operator/=(Numeric x)
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
inline Tensor3View& Tensor3View::operator+=(Numeric x)
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
inline Tensor3View& Tensor3View::operator-=(Numeric x)
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

/** Element-vise multiplication by another Tensor3. */
inline Tensor3View& Tensor3View::operator*=(const ConstTensor3View& x)
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

/** Element-vise division by another Tensor3. */
inline Tensor3View& Tensor3View::operator/=(const ConstTensor3View& x)
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

/** Element-vise addition of another Tensor3. */
inline Tensor3View& Tensor3View::operator+=(const ConstTensor3View& x)
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

/** Element-vise subtraction of another Tensor3. */
inline Tensor3View& Tensor3View::operator-=(const ConstTensor3View& x)
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

/** Element-vise multiplication by a Vector (acting like a 1-column Tensor3). */
inline Tensor3View& Tensor3View::operator*=(const ConstVectorView& x)
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

/** Element-vise division by a Vector (acting like a 1-column Tensor3). */
inline Tensor3View& Tensor3View::operator/=(const ConstVectorView& x)
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

/** Element-vise addition of a Vector (acting like a 1-column Tensor3). */
inline Tensor3View& Tensor3View::operator+=(const ConstVectorView& x)
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

/** Element-vise subtraction of a Vector (acting like a 1-column Tensor3). */
inline Tensor3View& Tensor3View::operator-=(const ConstVectorView& x)
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
			    const Range& rr,
			    const Range& cr) :
  ConstTensor3View(data, rr, cr)
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
			    const Range& pr, const Range& pc,
			    const Range& nr, const Range& nc) :
  ConstTensor3View(data,pr,pc,nr,nc)
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
inline void copy(ConstIterator2D origin,
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
inline void copy(Numeric x,
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

/** Constructor setting size. This constructor has to set the stride
    in the row range correctly! */
inline Tensor3::Tensor3(Index r, Index c) :
  Tensor3View( new Numeric[r*c],
	     Range(0,r,c),
	     Range(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Tensor3::Tensor3(Index r, Index c, Numeric fill) :
  Tensor3View( new Numeric[r*c],
	      Range(0,r,c),
	      Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x=mdata; x<mdata+r*c; ++x )
    *x = fill;
}

/** Copy constructor from Tensor3View. This automatically sets the size
    and copies the data. */
inline Tensor3::Tensor3(const ConstTensor3View& m) :
  Tensor3View( new Numeric[m.nrows()*m.ncols()],
	     Range( 0, m.nrows(), m.ncols() ),
	     Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Tensor3. This automatically sets the size
    and copies the data. */
inline Tensor3::Tensor3(const Tensor3& m) :
  Tensor3View( new Numeric[m.nrows()*m.ncols()],
	     Range( 0, m.nrows(), m.ncols() ),
	     Range( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty matrix, then it will have
  // 0 colunns. But this is used to initialize the stride of the row
  // Range! Thus, this method has to be consistent with the behaviour
  // of Range::Range. For now, Range::Range allows also stride 0.
  copy(m.begin(),m.end(),begin());
}

/** Assignment operator from another matrix. It is important that this
    operator exists. Otherwise the = operator seems to copy references
    instead of content in some cases. 

    The Behavior of this one is a bit special: If the size of the
    target Tensor3 is 0 then it will be automatically resized to match
    (this is needed to have the correct initialization for constructed
    classes that use the assignment operator to initialize their
    data). 
*/
inline Tensor3& Tensor3::operator=(const Tensor3& m)
{
  //  cout << "Tensor3 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "             n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid matrix, so we just
  // have to check one.
  if ( 0 == mrr.mextent )
    {
      // Adjust if previously empty.
      resize( m.mrr.mextent, m.mcr.mextent ); 
    }
  else
    {
      // Check that sizes are compatible:
      assert( mrr.mextent==m.mrr.mextent );
      assert( mcr.mextent==m.mcr.mextent );
    }

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment operator from scalar. Assignment operators also seem to
    be not inherited. */
inline Tensor3& Tensor3::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Assignment from a vector. This copies the data from a VectorView
    to this Tensor3View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor3View by
    setting its range. */
inline Tensor3& Tensor3::operator=(const ConstVectorView& v)
{
  // Check that sizes are compatible:
  assert( mrr.mextent==v.nelem() );
  assert( mcr.mextent==1         );
  //  dummy = ConstTensor3View(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstTensor3View dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Tensor3 is not
    initialized, so it will contain random values.*/
inline void Tensor3::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );

  if ( mrr.mextent!=r || mcr.mextent!=c )
    {
      delete mdata;
      mdata = new Numeric[r*c];

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


// Some general Matrix Vector functions:

/** Scalar product. The two vectors may be identical. */
inline Numeric operator*(const ConstVectorView& a, const ConstVectorView& b)
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

/** Matrix Vector multiplication. y = M*x. Note that the order is different
    from MTL, output comes first! Dimensions of y, M, and x must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping MatrixViews belonging
    to the same Matrix will lead to unpredictable results. */
inline void mult( VectorView y,
		  const ConstMatrixView& M,
		  const ConstVectorView& x )
{
  // Check dimensions:
  assert( y.nelem() == M.nrows() );
  assert( M.ncols() == x.nelem() );

  // Let's get a 1D iterator for y:
  const Iterator1D ye = y.end();
  Iterator1D       yi = y.begin();

  // ... a 2D iterator for M:
  ConstIterator2D  mi = M.begin();

  // ... and 1D iterators pointing to the start and end of x:
  const ConstIterator1D xs = x.begin();
  const ConstIterator1D xe = x.end();

  // This walks through the rows of y and M:
  for ( ; yi!=ye ; ++yi, ++mi )
    {
      // Compute the scalar product between this row of M and x: 
      Numeric dummy = 0;
      
      // Iterators for row of M:
      ConstIterator1D       ri = mi->begin();

      // ... and for x (always initialized to the start of x):
      ConstIterator1D       xi = xs;

      for ( ; xi!=xe ; ++xi, ++ri )
	{
	  dummy += (*ri) * (*xi);
	}

      *yi = dummy;
    }
}

/** Matrix multiplication. A = B*C. Note that the order is different
    from MTL, output comes first! Dimensions of A, B, and C must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping MatrixViews belonging
    to the same Matrix will lead to unpredictable results. */
inline void mult( MatrixView A,
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

/** Const version of transpose. */
inline ConstMatrixView transpose(ConstMatrixView m)
{
  return ConstMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
    transpose. The original is not changed! */
inline MatrixView transpose(MatrixView m)
{
  return MatrixView(m.mdata, m.mcr, m.mrr);
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

    \retval   y   the results of the function acting on each element of x
    \param    my_func a function (e.g., sqrt)
    \param    x   a vector */
inline void transform( VectorView y,
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

   \retval   y   the results of the function acting on each element of x
   \param    my_func a function (e.g., sqrt)
   \param    x   a matrix */
inline void transform( MatrixView y,
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
inline Numeric max(const ConstVectorView& x)
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
inline Numeric max(const ConstMatrixView& x)
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
inline Numeric min(const ConstVectorView& x)
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
inline Numeric min(const ConstMatrixView& x)
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


#endif    // matpackIII_h

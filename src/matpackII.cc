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

#include <vector>
#include "matpackI.h"
#include "matpackII.h"

/*
  Notes:

  There are two different ways to index:
  S.rw(3,4) = 1;                // Read and write
  cout << S.ro(3,4);            // Read only

  This distinction is necessary, because rw() creates elements if they
  don't already exist.

 */

// Functions for SparseView:
// -------------------------------

/** Returns the number of rows. */
Index SparseView::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
Index SparseView::ncols() const
{
  return mcr.mextent;
}

/** Returns the number of nonzero elements. */
Index SparseView::nnz() const
{
  return *(mcolptr->end()-1);
}

/** Plain index operator. This has to correctly handle two cases:

    1. The data element exists. In this case the operator acts similar
       to the const index operator in that the element is returned.

    2. The element does not exist. In this case it is created.
*/
Numeric& SparseView::rw(Index r, Index c)
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  // Convert to true indices. We use r and c, which are local
  // variables because we used call be value.
  r = mrr.mstart + r*mrr.mstride;
  c = mcr.mstart + c*mcr.mstride;

  // Get index of first data element of this column:
  Index i = (*mcolptr)[c];

  // Get index of first data element of next column. (This always works,
  // because mcolptr has one extra element pointing behind the last
  // column.)
  const Index end = (*mcolptr)[c+1];

  // See if we find an element with the right row index in this
  // range. We assume that the elements are sorted by ascending row
  // index. We use the index i as the loop counter, which we have
  // initialized above.
  for ( ; i<end; ++i )
    {
      Index rowi = (*mrowind)[i];
      if ( r <  rowi )
        break;
      if ( r == rowi )
        return (*mdata)[i];
    }

  // If we are here, then the requested data element does not yet
  // exist. 

  // We have to adjust the array of column pointers. The values
  // in all columns above the current one have to be increased by 1.
  for ( std::vector<Index>::iterator j = mcolptr->begin() + c + 1;
        j < mcolptr->end();
        ++j )
    ++(*j);

  // We have to insert the new element in *mrowind and *mdata just one
  // position before the index i. We can use std::insert to achieve
  // this. Because they return an iterator to the newly inserted
  // element, we can return directly from the second call.
  mrowind->insert( mrowind->begin()+i, r );
  return *( mdata->insert( mdata->begin()+i, 0 ) );
}

/** Plain const index operator. */
Numeric SparseView::ro(Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  // Convert to true indices:
  r = mrr.mstart + r*mrr.mstride;
  c = mcr.mstart + c*mcr.mstride;

  // Get index of first data element of this column:
  Index begin = (*mcolptr)[c];

  // Get index of first data element of next column. (This always works,
  // because mcolptr has one extra element pointing behind the last
  // column.)
  const Index end = (*mcolptr)[c+1];

  // See if we find an element with the right row index in this
  // range. We assume that the elements are sorted by ascending row
  // index. 
  for ( Index i=begin; i<end; ++i )
    {
      Index rowi = (*mrowind)[i];
      if ( r <  rowi )
        return 0;
      if ( r == rowi )
        return (*mdata)[i];
    }
  return 0;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Sparse. */
SparseView::SparseView() :
  mdata(NULL),
  mrowind(NULL),
  mcolptr(NULL),
  mrr(0,0),
  mcr(0,0)
{
  // Nothing to do here
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own SubMatrix part. The row range rr must have a
    stride to account for the length of one row. */
SparseView::SparseView(std::vector<Numeric> *data,
                                        std::vector<Index>   *rowind, 
                                        std::vector<Index>   *colptr,
                                        const Range&         rr,      
                                        const Range&         cr       ) :
  mdata(data),
  mrowind(rowind),
  mcolptr(colptr),
  mrr(rr),
  mcr(cr)
{
  // Nothing to do here.
}



// Functions for Sparse:
// ---------------------------

/** Default constructor. */
Sparse::Sparse() 
{
  // Nothing to do here
}

/** Constructor setting size.

    Elements *mdata and *mrowind have to grow later on, when we add
    data element. But *mcolptr always has the dimenson of the number
    of columns of the matrix plus one, so it is allocated
    directly. 

    Why is there an extra element in *mcolptr? We store also the index
    *behind* the last element of the last column. Or in other words
    the starting index that the next column *would* have. This just
    safes a litle time when computing indices. Also, this corresponds
    to the number of nonzero elements.
*/
Sparse::Sparse(Index r, Index c) :
  SparseView( new std::vector<Numeric>,
                   new std::vector<Index>,
                   new std::vector<Index>(c+1,0),
                   Range(0,r),
                   Range(0,c) )
{
  // Nothing to do here.
}

/** Copy constructor from another Sparse. This automatically
    sets the size and copies the data. */
Sparse::Sparse(const Sparse& m) :
  SparseView( new std::vector<Numeric>(*m.mdata),
                   new std::vector<Index>(*m.mrowind),
                   new std::vector<Index>(*m.mcolptr),
                   m.mrr,
                   m.mcr )
{
  // Nothing to do here. 
}


/** Destructor for Sparse. This is important, since Sparse
    uses new to allocate storage. */
Sparse::~Sparse()
{
  delete mdata;
  delete mrowind;
  delete mcolptr;
}

// Output operator for SparseView:

std::ostream& operator<<(std::ostream& os, const SparseView& v)
{
  for (size_t c=0; c<v.mcolptr->size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = (*v.mcolptr)[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = (*v.mcolptr)[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i )
        {
          Index r = (*v.mrowind)[i];

          // Now we know the true row r and column c. Convert to aparent r
          // and c using the Ranges mrr and mcr,

          Index ra = (r-v.mrr.mstart)/v.mrr.mstride;
          Index ca = (c-v.mcr.mstart)/v.mcr.mstride;

          // Are the aparent row ra and column ca inside the active range?
          if ( 0 <= ra &&
               ra < v.mrr.mextent &&
               0 <= ca &&
               ca < v.mcr.mextent     )
            {
              // Yes, they are! Let's output this element.
              os << setw(3) << ra << " "
                 << setw(3) << ca << " "
                 << setw(3) << (*v.mdata)[i] << "\n";
            }
        }
    }

  return os;
}

// General matrix functions

/** Matrix Vector multiplication. y = M*x. Note that the order is different
    from MTL, output comes first! Dimensions of y, M, and x must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping MatrixViews belonging
    to the same Matrix will lead to unpredictable results. In
    particular, this means that A and B must not be the same matrix! */
void mult( VectorView y,
           const SparseView& M,
           const ConstVectorView& x )
{
  // Check dimensions:
  assert( y.nelem() == M.nrows() );
  assert( M.ncols() == x.nelem() );

  // FIXME: Maybe this should be done with iterators as for Matrix
  // but then we need to define the 2D iterator for Sparse
  
  // Looping through every element of M, multiplying it with the
  // appropriate elements of x.
  for (size_t c=0; c<M.mcolptr->size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = (*M.mcolptr)[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = (*M.mcolptr)[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i )
        {
          Index r = (*M.mrowind)[i];

          // Now we know the true row r and column c. Convert to apparent r
          // and c using the Ranges mrr and mcr,

          Index ra = (r-M.mrr.mstart)/M.mrr.mstride;
          Index ca = (c-M.mcr.mstart)/M.mcr.mstride;
          
          // Convert apparent indices to true for the VectorViews
          //Index ry = ra*y.mrange.mstride + y.mrange.mstart;
          //Index cx = ca*x.mrange.mstride + x.mrange.mstart;

          // Are the aparent row ra and column ca inside the active range?
          if ( 0 <= ra &&
               ra < M.mrr.mextent &&
               0 <= ca &&
               ca < M.mcr.mextent     )
            {
              // Yes, they are! Let's compute this element.
              // y[i] = M(i,j) * x[j]
              //y[ry] = (*M.mdata)[i] * x[cx];
              y[ra] = (*M.mdata)[i] * x[ca];

              // FIXME: apparent indices or true (in the vectors)
              /*
              os << setw(3) << ra << " "
                 << setw(3) << ca << " "
                 << setw(3) << (*v.mdata)[i] << "\n";
              */
            }
        }
    }
}

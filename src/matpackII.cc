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
#include <algorithm>
#include <set>
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

Numeric& SparseView::rw (Index r, Index c)
{ return (*this)(r, c); }

/** Plain index operator. This has to correctly handle two cases:

    1. The data element exists. In this case the operator acts similar
       to the const index operator in that the element is returned.

    2. The element does not exist. In this case it is created.
*/
Numeric& SparseView::operator ()(Index r, Index c)
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

Numeric SparseView::ro (Index r, Index c) const
{ return (*this)(r, c); }

/** Plain const index operator. */
Numeric SparseView::operator() (Index r, Index c) const
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
    data element. But *mcolptr always has the dimension of the number
    of columns of the matrix plus one, so it is allocated
    directly. 

    Why is there an extra element in *mcolptr? We store also the index
    *behind* the last element of the last column. Or in other words
    the starting index that the next column *would* have. This just
    safes a little time when computing indices. Also, this corresponds
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

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Sparse is not
    initialised so it will be empty.*/
void Sparse::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );
  if ( mrr.mextent!=r || mcr.mextent!=c )
    {
      delete mdata;
      mdata = new std::vector<Numeric>;
      delete mrowind;
      mrowind = new std::vector<Index>;
      delete mcolptr;
      mcolptr = new std::vector<Index>(c+1,0);

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
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

          // Now we know the true row r and column c. Convert to apparent r
          // and c using the Ranges mrr and mcr,

          Index ra = (r-v.mrr.mstart)/v.mrr.mstride;
          Index ca = (c-v.mcr.mstart)/v.mcr.mstride;

          // Are the apparent row ra and column ca inside the active range?
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

          // Are the apparent row ra and column ca inside the active range?
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

/** SparseMatrix - Matrix multiplication. A = B*C, where B is sparse.
    Note that the order is different from MTL, output comes first!
    Dimensions of A, B, and C must match. No memory reallocation takes
    place, only the data is copied. Using this function on overlapping
    MatrixViews belonging to the same Matrix will lead to unpredictable
    results. In particular, this means that A and C must not be the
    same matrix! */
void mult( MatrixView A,
           const SparseView B,
           const MatrixView C )
{
  // Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Set all elements of A to zero
  A = 0.0;

  // Loop through the element of B
  for (size_t c=0; c<B.mcolptr->size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = (*B.mcolptr)[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = (*B.mcolptr)[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i ) {
        Index r = (*B.mrowind)[i];

        // Now we know the true row r and column c. Convert to apparent r
        // and c using the Ranges mrr and mcr,
        Index ra = (r-B.mrr.mstart)/B.mrr.mstride;
        Index ca = (c-B.mcr.mstart)/B.mcr.mstride;

        // Are the apparent row ra and column ca inside the active range?
        if ( 0 <= ra &&
             ra < B.mrr.mextent &&
             0 <= ca &&
             ca < B.mcr.mextent     )
        {
          // Yes, they are! Multiply it with corresponding column in C
          // and add the product to the right element in A
          for (Index j=0; j<C.ncols(); j++) {
            A(ra,j) += B.ro(ra,ca) * C(ca,j);
          }
        }
      }
    }
}

/** Transpose of sparse matrix

    2003-04-04  Mattias Ekström
*/
void transpose( SparseView A,
                const SparseView B )
{
  // Check that sizes match
  assert( A.nrows() == B.ncols() );
  assert( A.ncols() == B.nrows() );

  // Create a vector with all row indices, sorted in strict ascending order
  std::vector<Index> rowind = *B.mrowind;
  sort(rowind.begin(), rowind.end());
  std::vector<Index>::iterator last = unique(rowind.begin(), rowind.end());
  rowind.erase(last, rowind.end());

  // Loop through the existing rows of B and check them vs columns of C
  for (Index l=0; l<rowind.size(); l++) {
    Index i = rowind[l];
    
    // Loop through columns and get the values for the specific row
    for (Index j=0; j<B.ncols(); j++) {
      // Get index of first data element of this column:
      Index begin = (*B.mcolptr)[j];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      Index end = (*B.mcolptr)[j+1];

      // If row index is within the span of this column, search for it
      if (i>=(*B.mrowind)[begin] && i<=(*B.mrowind)[end-1]) {
        for (Index k=begin; k<end; ++k) {
          if ( i == (*B.mrowind)[k] )
            A.rw(j,i) = B.ro(i,j);
        }
      }
    }
  }
  // Should rowind be destructed?
}


/** Sparse - Sparse multiplication. A = B*C, where result A is sparse.

    Note that the order is different from MTL, output comes first!
    Dimensions of A, B, and C must match. No memory reallocation takes
    place, only the data is copied. Using this function on overlapping
    SparseViews belonging to the same Sparse will lead to unpredictable
    results. In particular, this means that A and C or A and B must not
    be the same matrix!

    2003-04-04  Mattias Ekström
*/
void mult( SparseView A,
           const SparseView B,
           const SparseView C )
{
  //Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  //Transpose B to simplify multiplication algorithm
  Sparse Bt(B.ncols(), B.nrows());
  transpose(Bt,B);

  //Loop over columns in C and multiply them with every column in Bt
  for (size_t cC=0; cC<C.mcolptr->size()-1; ++cC) {
    //Get row indices of this column
    Index beginC = (*C.mcolptr)[cC];
    Index endC = (*C.mcolptr)[cC+1];

    for (size_t cBt=0; cBt<Bt.mcolptr->size()-1; ++cBt) {
      //Get row indices for this column too
      Index beginBt = (*Bt.mcolptr)[cBt];
      Index endBt = (*Bt.mcolptr)[cBt+1];

      /*
      cout << "Test: "<<(*Bt.mcolptr)[cBt]<<":"<<(*Bt.mcolptr)[cBt+1]<<"-"<<(*C.mcolptr)[cC]<<":"<<(*C.mcolptr)[cC+1]<<" ";
      if( endBt-beginBt!=0 )
        cout << "Bt ok";
      if( endC-beginC!=0 )
        cout << ", C ok";
      if( (*Bt.mrowind)[endBt-1]>=(*C.mrowind)[beginC] )
        cout << ", Bt(end)>C(first) ok";
      if( (*Bt.mrowind)[beginBt]<=(*C.mrowind)[endC-1] )
        cout << ", Bt(first)<C(end) ok";
      cout << "\n";
      */

      //Check that the columns are non-empty, ...
      if ( endBt-beginBt!=0 && endC-beginC!=0
          // (NB: last index actually points to next columns first)
          // that they are overlapping and ...
          && (*Bt.mrowind)[endBt-1]>=(*C.mrowind)[beginC]
          && (*Bt.mrowind)[beginBt]<=(*C.mrowind)[endC-1]
          // that this special case of overlapping is not included.
          /*&& beginBt!=endC && beginC!=endBt*/ ) {
        Numeric tempA=0.0;

        //Go through columns and find matching indices
        Index i=beginBt, j=beginC;
        while ( j<endC && i<endBt ) {
          //cout<<"B("<<cBt<<","<<i<<")*C("<<j<<","<<cC<<")="<<Bt.ro(i,cBt)<<"*"<<C.ro(j,cC)<<"="<<Bt.ro(i,cBt)*C.ro(j,cC)<<endl;
          cout <<"i="<<(*Bt.mrowind)[i]<<",j="<<(*C.mrowind)[j]<<",";
          if ((*C.mrowind)[j]>(*Bt.mrowind)[i]) {
            i++;
            cout << "i++,";
          } else if ((*C.mrowind)[j]<(*Bt.mrowind)[i]) {
            j++;
            cout << "j++,";
          } else {
            tempA += (*Bt.mdata)[i] * (*C.mdata)[j];
            i++;
            j++;
          }
        }
        cout << " tempA " << tempA << "\n";

        //Did we get a sum?
        if (tempA!=0.0) {
          //Yes, write it to product A
          A.rw(cBt,cC) = tempA;
        }
      }
    }
  }
}

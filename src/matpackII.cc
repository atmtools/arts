/* Copyright (C) 2001, 2002, 2003
   Stefan Buehler <sbuehler@uni-bremen.de>
   Mattias Ekstroem <ekstrom@rss.chalmers.se>

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
  \file   matpackII.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003
  
  \brief  Implementation of sparse matrices.
  
  Notes:

  There are two different ways to index: 
  S.rw(3,4) = 1;                // Read and write
  cout << S.ro(3,4);            // Read only

  This distinction is necessary, because rw() creates elements if they
  don't already exist.

  The normal index operator "()" correspondes to ro, so "S(3,4)" is
  the same as S.ro(3,4).
*/

// #include <vector>
// #include <iostream>
// #include <algorithm>
// #include <set>
#include "matpackII.h"


// Simple member Functions
// ----------------

//! Returns the number of rows. 
Index Sparse::nrows() const
{
  return mrr;
}

//! Returns the number of columns. 
Index Sparse::ncols() const
{
  return mcr;
}

//! Returns the number of nonzero elements. 
Index Sparse::nnz() const
{
  return *(mcolptr->end()-1);
}

// Index Operators
// ---------------

//! Read and write index operator.
/*! 
  This has to correctly handle two cases:

  1. The data element exists. In this case the operator acts similar
     to the const index operator in that the element is returned.
     
  2. The element does not exist. In this case it is created.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.
  
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
Numeric& Sparse::rw(Index r, Index c)
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr );
  assert( c<mcr );

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

//! Plain index operator.
/*! 
  This is the same as the .ro index operator.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.
  
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
Numeric Sparse::operator() (Index r, Index c) const
{ return this->ro(r, c); }

//! Read only index operator.
/*! 
  This has to correctly handle two cases:

  1. The data element exists. In this case the element is returned.
     
  2. The element does not exist. In this case the value 0 is returned.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices, or zero.
  
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 */
Numeric Sparse::ro (Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr );
  assert( c<mcr );

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


// Constructors
// ------------

//! Default constructor.
/*!
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
Sparse::Sparse() :
  mdata(NULL),
  mrowind(NULL),
  mcolptr(NULL),
  mrr(0),
  mcr(0)
{
  // Nothing to do here
}

//! Constructor setting size.
/*! 
  Elements *mdata and *mrowind have to grow later on, when we add
  data elements. But *mcolptr always has the dimension of the number
  of columns of the matrix plus one, so it is allocated
  directly. (And properly initialized to zero.)

  Why is there an extra element in *mcolptr? We store also the index
  *behind* the last element of the last column. Or in other words
  the starting index that the next column *would* have. This just
  safes a little time when computing indices. Also, this corresponds
  to the number of nonzero elements.
  
  \param r Row dimension of new sparse matrix.
  \param c Column dimension of new sparse matrix.

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
Sparse::Sparse(Index r, Index c) :
  mdata(new std::vector<Numeric>),
  mrowind(new std::vector<Index>),
  mcolptr(new std::vector<Index>(c+1,0)),    
  mrr(r),
  mcr(c)
{
  // Nothing to do here.
}


//! Copy constructor from another Sparse.
/*! 
  This automatically sets the size and copies the data.
  
  \param m The other Sparse to copy from.

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
Sparse::Sparse(const Sparse& m) :
  mdata(new std::vector<Numeric>(*m.mdata)),
  mrowind(new std::vector<Index>(*m.mrowind)),
  mcolptr(new std::vector<Index>(*m.mcolptr)),    
  mrr(m.mrr),
  mcr(m.mcr)
{
  // Nothing to do here. 
}


//! Destructor for Sparse.
/*! 
  This is important, since Sparse uses new to allocate storage.

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
Sparse::~Sparse()
{
  delete mdata;
  delete mrowind;
  delete mcolptr;
}

//! Resize function.
/*!
  If the size is already correct this function does nothing. 

  All data is lost after resizing! The new Sparse is not initialised
  so it will be empty.

  \param r New row dimension.
  \param c New column dimension.

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
void Sparse::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );
  if ( mrr!=r || mcr!=c )
    {
      delete mdata;
      mdata = new std::vector<Numeric>;
      delete mrowind;
      mrowind = new std::vector<Index>;
      delete mcolptr;
      mcolptr = new std::vector<Index>(c+1,0);

      mrr = r;
      mcr = c;
    }
}

//! Output operator for Sparse.
/*!   
  \param os Output stream.
  \param v Sparse matrix to print.

  \return Output stream.

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
std::ostream& operator<<(std::ostream& os, const Sparse& v)
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
          // Get row index:
          Index r = (*v.mrowind)[i];

          // Output everything:
          os << setw(3) << r << " "
             << setw(3) << c << " "
             << setw(3) << (*v.mdata)[i] << "\n";
        }
    }

  return os;
}

// General matrix functions

//! Sparse matrix - Vector multiplication.
/*!
  This calculates the product

  y = M*x, where M is sparse.

  Output comes first! 

  Dimensions of y, M, and x must match. No memory reallocation takes
  place, only the data is copied.

  \param y Output: The multiplication result.
  \param M Matrix for multiplication (sparse).
  \param x Vector for multiplication.

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
void mult( VectorView y,
           const Sparse& M,
           ConstVectorView x )
{
  // Check dimensions:
  assert( y.nelem() == M.nrows() );
  assert( M.ncols() == x.nelem() );

  // Initialize y to all zeros:
  y = 0.0;

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
          // Get row index:
          Index r = (*M.mrowind)[i];

          // Compute this element:
          y[r] += (*M.mdata)[i] * x[c];
        }
    }
}

//! SparseMatrix - Matrix multiplication.
/*!
  Calculates the matrix product:

  A = B*C, where B is sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.
 
  \param A Output: Result matrix (full).
  \param B First matrix to multiply (sparse).
  \param C Second matrix to multiply (full).

  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Tue Jul 15 15:05:40 2003 
*/
void mult( MatrixView A,
           const Sparse& B,
           ConstMatrixView C )
{
  // Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Set all elements of A to zero:
  A = 0.0;

  // Loop through the elements of B:
  for (size_t c=0; c<B.mcolptr->size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = (*B.mcolptr)[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = (*B.mcolptr)[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i )
        {
          // Get row index:
          Index r = (*B.mrowind)[i];

          // Multiply this element with the corresponding row of C and
          // add the product to the right row of A
          for (Index j=0; j<C.ncols(); j++)
            {
              /* Conceptually:
                 A(r,j) += B.ro(r,c) * C(c,j);
                 But we don't need to use the index operator, because
                 we have the right element right here: */
              A(r,j) += (*B.mdata)[i] * C(c,j);
            }
        }
    }
}


//! Transpose of sparse matrix
/*!
  \param A Output: Transposed matrix.
  \param B Original matrix.

  \author Mattias Ekstroem
  \date   2003-04-04  
*/
void transpose( Sparse& A,
                const Sparse& B )
{
  /*
    FIXME: Mattias, since this is a friend function, you can access the
    internal data elements directly. Then you don't have to use the
    .rw() operator to insert elements in the new matrix, but can reserve
    the right amount of space right away. 

    I don't understand your code, so I'm not attempting to fix it. I've
    just commented it out for now. Just use my other functions as
    samples and formulate your algorithm.

    I strongly recommend to use more explicit comments. ;-)

    - Stefan
  */

  exit(1);

  /*
  // Check that sizes match
  assert( A.nrows() == B.ncols() );
  assert( A.ncols() == B.nrows() );

  // Create a vector with all row indices, sorted in strict ascending order
  std::vector<Index> rowind = *B.mrowind;
  sort(rowind.begin(), rowind.end());
  std::vector<Index>::iterator last = unique(rowind.begin(), rowind.end());
  rowind.erase(last, rowind.end());

  // Loop through the existing rows of B and check them vs columns of C
  for (Index l=0; l<(signed)rowind.size(); l++) {
    Index i = (rowind[l]-B.mrr.get_start ())/B.mrr.get_stride ();

    // Loop through columns and get the values for the specific row
    for (Index j=0; j<B.ncols(); j++) {
      // Get index of first data element of this column:
      Index begin = (*B.mcolptr)[j];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      Index end = (*B.mcolptr)[j+1];

      // If row index is within the span of this column, search for it
      Index firstB = ((*B.mrowind)[begin]-B.mrr.get_start ())/B.mrr.get_stride ();
      Index lastB = ((*B.mrowind)[end-1]-B.mrr.get_start ())/B.mrr.get_stride ();
      if (i>=firstB && i<=lastB) {
        for (Index k=begin; k<end; ++k) {
          if ( i == ((*B.mrowind)[k]-B.mrr.get_start ())/B.mrr.get_stride ())
            A.rw(j,i) = B.ro(i,j);
        }
      }
    }
  }
  */
}

//! Sparse - Sparse multiplication.
/*!
  Calculates A = B*C, where result A is sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.
  
  \param A Output: Result matrix.
  \param B First product matrix.
  \param C Second product matrix.

  \author Mattias Ekstroem
  \date   2003-04-04  
*/
void mult( Sparse& A,
           const Sparse& B,
           const Sparse& C )
{
  /*
    FIXME: Mattias, as for transpose, I don't understand your code, so
    I'm commenting it out for now. ;-)

    - Stefan
  */

  exit(1);

  /*
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

//       cout << "Test: "<<(*Bt.mcolptr)[cBt]<<":"<<(*Bt.mcolptr)[cBt+1]<<"-"<<(*C.mcolptr)[cC]<<":"<<(*C.mcolptr)[cC+1]<<" ";
//       if( endBt-beginBt!=0 )
//         cout << "Bt ok";
//       if( endC-beginC!=0 )
//         cout << ", C ok";
//       if( (*Bt.mrowind)[endBt-1]>=(*C.mrowind)[beginC] )
//         cout << ", Bt(end)>C(first) ok";
//       if( (*Bt.mrowind)[beginBt]<=(*C.mrowind)[endC-1] )
//         cout << ", Bt(first)<C(end) ok";
//       cout << "\n";

      //Check that the columns are non-empty, ...
      Index firstBt = ((*Bt.mrowind)[beginBt]-Bt.mrr.get_start ())
        / Bt.mrr.get_stride ();
      Index lastBt = ((*Bt.mrowind)[endBt-1]-Bt.mrr.get_start ())
        / Bt.mrr.get_stride ();
      Index firstC = ((*C.mrowind)[beginC]-C.mrr.get_start ())
        / C.mrr.get_stride ();
      Index lastC = ((*C.mrowind)[endC-1]-C.mrr.get_start ())
        / C.mrr.get_stride ();
      if ( endBt-beginBt!=0 && endC-beginC!=0
          // (NB: last index actually points to next columns first)
          // that they are overlapping and ...
          && lastBt>=firstC && firstBt<=lastC
          //&& (*Bt.mrowind)[endBt-1]>=(*C.mrowind)[beginC]
          //&& (*Bt.mrowind)[beginBt]<=(*C.mrowind)[endC-1]
          // that this special case of overlapping is not included.
          // && beginBt!=endC && beginC!=endBt ) {
        Numeric tempA=0.0;

        //Go through columns and find matching indices
        Index i=beginBt, j=beginC;
        while ( j<endC && i<endBt ) {
          //cout<<"B("<<cBt<<","<<i<<")*C("<<j<<","<<cC<<")="<<Bt.ro(i,cBt)<<"*"<<C.ro(j,cC)<<"="<<Bt.ro(i,cBt)*C.ro(j,cC)<<endl;
          //cout <<"i="<<(*Bt.mrowind)[i]<<",j="<<(*C.mrowind)[j]<<",";
          Index c = ((*C.mrowind)[j]-C.mrr.get_start ())/C.mrr.get_stride ();
          Index bt =((*Bt.mrowind)[i]-Bt.mrr.get_start ())/Bt.mrr.get_stride ();
          //if ((*C.mrowind)[j]>(*Bt.mrowind)[i]) {
          if (c>bt) {
            i++;
            //cout << "i++,";
          //} else if ((*C.mrowind)[j]<(*Bt.mrowind)[i]) {
          } else if (c<bt) {
            j++;
            //cout << "j++,";
          } else {
            tempA += (*Bt.mdata)[i] * (*C.mdata)[j];
            i++;
            j++;
          }
        }
        //cout << " tempA " << tempA << "\n";

        //Did we get a sum?
        if (tempA!=0.0) {
          //Yes, write it to product A
          A.rw(cBt,cC) = tempA;
        }
      }
    }
  }
*/
}

//! Sparse - Sparse multiplication, version 2.
/*!
  This is a second version of the sparse-sparse multiplication
  algorithm. The aim is to hold down the memory usage created when
  making a transposed copy of B.
  
  \param A Output: Result matrix.
  \param B First product matrix.
  \param C Second product matrix.

  \author Mattias Ekstroem
  \date   2003-06-27  
*/
void mult2( Sparse A,
           const Sparse B,
           const Sparse C )
{
  exit(1);

  /*
  //Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Create a vector with all row indices, sorted in strict ascending order
  std::vector<Index> rowind = *B.mrowind;
  sort(rowind.begin(), rowind.end());
  std::vector<Index>::iterator last = unique(rowind.begin(), rowind.end());
  rowind.erase(last, rowind.end());

  // Loop through the existing rows of B and check them vs columns of C
  for (Index l=0; l<(signed)rowind.size(); l++) {
    Index i = (rowind[l]-B.mrr.mstart)/B.mrr.mstride;

    // Loop through columns and get the values for the specific row
    for (Index j=0; j<B.ncols(); j++) {
      // Get index of first data element of this column:
      Index begin = (*B.mcolptr)[j];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      Index end = (*B.mcolptr)[j+1];

      // If row index is within the span of this column, search for it
      Index firstB = ((*B.mrowind)[begin]-B.mrr.mstart)/B.mrr.mstride;
      Index lastB = ((*B.mrowind)[end-1]-B.mrr.mstart)/B.mrr.mstride;
      if (i>=firstB && i<=lastB) {
        for (Index k=begin; k<end; ++k) {
          if ( i == ((*B.mrowind)[k]-B.mrr.mstart)/B.mrr.mstride )
            A.rw(j,i) = B.ro(i,j);
        }
      }
    }
  }
  // Should rowind be destructed?
*/
}

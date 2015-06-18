/* Copyright (C) 2001-2012
   Stefan Buehler <sbuehler@ltu.se>
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
  \author Stefan Buehler <sbuehler@ltu.se>
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
#include <algorithm>
#include <set>
#include <iostream>             // For debugging.
#include <cmath>
#include <iterator>
#include "matpackII.h"

using std::vector;
using std::setw;


// Simple member Functions
// ----------------

//! Returns true if variable size is zero.
bool Sparse::empty() const
{
    return (nrows() == 0 || ncols() == 0);
}

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
  return *(mcolptr.end()-1);
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
  
  \author Stefan Buehler <sbuehler@ltu.se>
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
  Index i = mcolptr[c];

  // Get index of first data element of next column. (This always works,
  // because mcolptr has one extra element pointing behind the last
  // column.)
  const Index end = mcolptr[c+1];

  // See if we find an element with the right row index in this
  // range. We assume that the elements are sorted by ascending row
  // index. We use the index i as the loop counter, which we have
  // initialized above.
  for ( ; i<end; ++i )
    {
      Index rowi = mrowind[i];
      if ( r <  rowi )
        break;
      if ( r == rowi )
        return mdata[i];
    }

  // If we are here, then the requested data element does not yet
  // exist. 

  // We have to adjust the array of column pointers. The values
  // in all columns above the current one have to be increased by 1.
  for ( vector<Index>::iterator j = mcolptr.begin() + c + 1;
        j < mcolptr.end();
        ++j )
    ++(*j);

  // We have to insert the new element in *mrowind and *mdata just one
  // position before the index i. We can use insert to achieve
  // this. Because they return an iterator to the newly inserted
  // element, we can return directly from the second call.
  mrowind.insert( mrowind.begin()+i, r );
  return *( mdata.insert( mdata.begin()+i, 0 ) );
}

//! Plain index operator.
/*! 
  This is the same as the .ro index operator.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.
  
  \author Stefan Buehler <sbuehler@ltu.se>
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
  
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003 */
Numeric Sparse::ro (Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr );
  assert( c<mcr );

  // Get index of first data element of this column:
  Index begin = mcolptr[c];

  // Get index of first data element of next column. (This always works,
  // because mcolptr has one extra element pointing behind the last
  // column.)
  const Index end = mcolptr[c+1];

  // See if we find an element with the right row index in this
  // range. We assume that the elements are sorted by ascending row
  // index. 
  for ( Index i=begin; i<end; ++i )
    {
      Index rowi = mrowind[i];
      if ( r <  rowi )
        return 0;
      if ( r == rowi )
        return mdata[i];
    }
  return 0;
}


// Constructors
// ------------

//! Default constructor.
/*!
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003 
*/
Sparse::Sparse() :
  mdata(),
  mrowind(),
  mcolptr(),
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

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003 
*/
Sparse::Sparse(Index r, Index c) :
  mdata(),
  mrowind(),
  mcolptr(c+1,0),
  mrr(r),
  mcr(c)
{
  // Nothing to do here.
}


//! Insert row function
/*!
  Inserts a Vector as row of elements at the given position.
  
  The row index must agree with the size of the matrix. This
  function can not be used to expand the matrix.
  Only non-zero values will be stored. If the destination row
  already exist it will be overwritten.

  \param r Where to insert the row
  \param v Vector to be inserted.

  \author Mattias Ekstr?m <ekstrom@rss.chalmers.se>
  \date   2003-08-11
*/
void Sparse::insert_row(Index r, Vector v)
{
  // Check if the row index and the Vector length are valid
  assert( 0<=r );
  assert( r<mrr );
  assert( v.nelem()==mcr );

  // Calculate number of non-zero elements in Vector v.
  Index vnnz=0;
  for (Index i=0; i<v.nelem(); i++)
    {
      if (v[i]!=0)
        {
          vnnz++;
        }
    }

  // Count number of already existing elements in this row. Create
  // reference mrowind and mdata vector that copies of the real mrowind and
  // mdata. Resize the real mrowind and mdata to the correct output size.
  Index rnnz = vnnz - count(mrowind.begin(),mrowind.end(),r);

  vector<Index> mrowind_ref(mrowind.size());
  copy(mrowind.begin(), mrowind.end(), mrowind_ref.begin());

  vector<Numeric> mdata_ref(mdata.size());
  copy(mdata.begin(), mdata.end(), mdata_ref.begin());

  mrowind.resize(mrowind_ref.size()+rnnz);
  mdata.resize(mdata_ref.size()+rnnz);

  // Create iterators to the output vectors to keep track of current
  // positions.
  vector<Index>::iterator mrowind_it = mrowind.begin();
  vector<Numeric>::iterator mdata_it = mdata.begin();

  // Create a variable to store the change to mcolptr for each run
  Index colptr_mod = 0;

  // Loop through Vector v and insert the non-zero elements.
  for (Index i=0; i<v.nelem(); i++)
    {
      // Get mdata- and mrowind iterators to start and end of this
      // (the i:th) reference column.
      vector<Numeric>::iterator dstart = mdata_ref.begin()+mcolptr[i];
      vector<Numeric>::iterator dend = mdata_ref.begin()+mcolptr[i+1];
      vector<Index>::iterator rstart = mrowind_ref.begin()+mcolptr[i];
      vector<Index>::iterator rend = mrowind_ref.begin()+mcolptr[i+1];

      // Apply mcolptr change now that we have the iterators to the
      // data and row indices.
      mcolptr[i] = colptr_mod;

      if (v[i]!=0)
        {
          // Check if r exist within this column, and get iterator for
          // mdata
          vector<Index>::iterator rpos = find(rstart,rend,r);
          vector<Numeric>::iterator dpos = dstart+(rpos-rstart);
          if (rpos!=rend)
            {
              // The index was found, replace the value in mdata with the
              // value from v.
              *dpos = v[i];

              // Copy this column to the ouput vectors.
              copy(rstart,rend,mrowind_it);
              copy(dstart,dend,mdata_it);

              // Adjust the position iterators accordingly.
              mrowind_it += rend-rstart;
              mdata_it += dend-dstart;

              // Set the mcolptr step, for next loop
              colptr_mod = mcolptr[i]+(rend-rstart);
            }
          else
            {
              // The row index was not found, look for the first index
              // greater than r
              rpos = find_if(rstart,rend,bind2nd(std::greater<Index>(),r));
              dpos = dstart+(rpos-rstart);

              // Copy the first part of the column to the output vector.
              copy(rstart,rpos,mrowind_it);
              copy(dstart,dpos,mdata_it);

              // Make sure mrowind_it and mdata_it points at the first
              // 'empty' position.
              mrowind_it += rpos-rstart;
              mdata_it += dpos-dstart;

              // Insert the new value from v in mdata and the row index r in
              // mrowind.
              *mrowind_it = r;
              *mdata_it = v[i];

              // Again, make sure mrowind_it and mdata_it points at the
              // first 'empty' position.
              mrowind_it++;
              mdata_it++;

               // Copy the rest of this column to the output vectors.
              copy(rpos,rend,mrowind_it);
              copy(dpos,dend,mdata_it);

              // Adjust the iterators a last time
              mrowind_it += rend-rpos;
              mdata_it += dend-dpos;

              // Set the mcolptr step, for next loop
              colptr_mod = mcolptr[i]+(rend-rstart+1);
            }
        }
      else
        {
          // Check if r exist within this column, and get iterator for
          // mdata
          vector<Index>::iterator rpos = find(rstart,rend,r);
          vector<Numeric>::iterator dpos = dstart+(rpos-rstart);
          if (rpos!=rend)
            {
              // The index was found and we use remove_copy to copy the rest
              // of the column to mrowind_new and mdata_new.
              remove_copy(rstart,rend,mrowind_it, *rpos);
              remove_copy(dstart,dend,mdata_it, *dpos);

              // Increase the mrowind_it and mdata_it
              mrowind_it += rend-rstart-1;
              mdata_it += dend-dstart-1;

              // Set the mcolptr step, for next loop
              colptr_mod = mcolptr[i]+(rend-rstart-1);
            }
          else
            {
              // The row index was not found, all we need is to copy this
              // columns to the output mrowind and mdata vectors.
              copy(rstart,rend,mrowind_it);
              copy(dstart,dend,mdata_it);

              // Adjust the mrowind_it and mdata_it iterators
              mrowind_it += rend-rstart;
              mdata_it += dend-dstart;

              // Set the mcolptr step, for next loop
              colptr_mod = mcolptr[i]+(rend-rstart);
            }
        }
    }
  // Apply mcolptr change for the one extra mcolptr element.
  *(mcolptr.end()-1) = colptr_mod;
}

//! Make Identity matrix
/*!
  This functions sets the Sparse matrix to be an identity matrix.

  The matrix will be remade to fit the given number of rows and columns.

  \param r New row dimension.
  \param c New column dimension.

  \author Mattias Ekstr?m
  \date   2003-12-05
*/
void Sparse::make_I( Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );

  // First get number of ones in the identity matrix
  Index n = std::min(r,c);

  // Remake and assign values to vectors
  mcolptr = vector<Index>( c+1, n);
  mrowind = vector<Index>(n);
  mdata = vector<Numeric>(n,1.0);

  // Loop over number of ones and assign values [0:n-1] to mcolptr and
  // mrowind
  vector<Index>::iterator rit = mrowind.begin();
  vector<Index>::iterator cit = mcolptr.begin();
  for( Index i=0; i<n; i++ ) {
    *rit++ = i;
    *cit++ = i;
  }

  mrr = r;
  mcr = c;
}

//! Resize function.
/*!
  If the size is already correct this function does nothing.

  All data is lost after resizing! The new Sparse is not initialised
  so it will be empty.

  \param r New row dimension.
  \param c New column dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
void Sparse::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );
  if ( mrr!=r || mcr!=c )
    {
      mdata.resize(0);
      mrowind.resize(0);
      mcolptr = vector<Index>(c+1,0);

      mrr = r;
      mcr = c;
    }
}


//! Output operator for Sparse.
/*!
  \param os Output stream.
  \param v Sparse matrix to print.

  \return Output stream.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
std::ostream& operator<<(std::ostream& os, const Sparse& v)
{
  for (size_t c=0; c<v.mcolptr.size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = v.mcolptr[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = v.mcolptr[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i )
        {
          // Get row index:
          Index r = v.mrowind[i];

          // Output everything:
          os << setw(3) << r << " "
             << setw(3) << c << " "
             << setw(3) << v.mdata[i] << "\n";
        }
    }

  return os;
}

// General matrix functions

//! Absolute value of sparse matrix elements
/*!
  Computes the absolute values of the elements in sparse matrix B.

  The output matrix A must have been initialized with the correct size.

  \param A Output: Absolute value matrix.
  \param B Original matrix.

  \author Mattias Ekstrom
  \date   2005-03-21
*/
void abs(       Sparse& A,
          const Sparse& B )
{
  // Check dimensions
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == B.ncols() );

  // Here we allocate memory for the A.mdata vector so that it matches the
  // input matrix, and then store the absolute values of B.mdata in it.
  A.mdata.resize( B.mdata.size() );
  Index end = B.mdata.size();
  for (Index i=0; i<end; i++)
    {
      A.mdata[i] = fabs(B.mdata[i]);
    }
  
  // The column pointer and row index vectors are copies of the input  
  A.mcolptr.resize( B.mcolptr.size() );
  copy(B.mcolptr.begin(), B.mcolptr.end(), A.mcolptr.begin());
  A.mrowind.resize( B.mrowind.size() );
  copy(B.mrowind.begin(), B.mrowind.end(), A.mrowind.begin());
}


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

  \author Stefan Buehler <sbuehler@ltu.se>
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
  for (size_t c=0; c<M.mcolptr.size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = M.mcolptr[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = M.mcolptr[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i )
        {
          // Get row index:
          Index r = M.mrowind[i];

          // Compute this element:
          y[r] += M.mdata[i] * x[c];
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

  \author Stefan Buehler <sbuehler@ltu.se>
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
  for (size_t c=0; c<B.mcolptr.size()-1; ++c)
    {
      // Get index of first data element of this column:
      Index begin = B.mcolptr[c];

      // Get index of first data element of next column. (This always works,
      // because mcolptr has one extra element pointing behind the last
      // column.)
      const Index end = B.mcolptr[c+1];

      // Loop through the elements in this column:
      for ( Index i=begin; i<end; ++i )
        {
          // Get row index:
          Index r = B.mrowind[i];

          // Multiply this element with the corresponding row of C and
          // add the product to the right row of A
          for (Index j=0; j<C.ncols(); j++)
            {
              /* Conceptually:
                 A(r,j) += B.ro(r,c) * C(c,j);
                 But we don't need to use the index operator, because
                 we have the right element right here: */
              A(r,j) += B.mdata[i] * C(c,j);
            }
        }
    }
}


//! Transpose of sparse matrix
/*!
  Computes the transpose of the sparse matrix B.

  The output matrix A must have been initialized with the correct size.

  \param A Output: Transposed matrix.
  \param B Original matrix.

  \author Mattias Ekstroem
  \date   2003-08-05
*/
void transpose( Sparse& A,
                const Sparse& B )
{
  // Check dimensions
  assert( A.nrows() == B.ncols() );
  assert( A.ncols() == B.nrows() );
  if ( B.mdata.size() == 0) return;

  // Here we allocate memory for the A.mdata and A.mrowind vectors so that
  // it matches the input matrix. (NB: that mdata and mrowind already has
  // the same size.)
  A.mdata.resize( B.mdata.size() );
  A.mrowind.resize( B.mrowind.size() );

  // Find the minimum and maximum existing row number in B.
  // (This maybe unnecessary in our case since we mostly treat diagonally
  // banded matrices where the minimum and maximum existing row numbers
  // coincide with the border numbers of the matrix.)
  vector<Index>::const_iterator startrow, stoprow;
  startrow = min_element( B.mrowind.begin(), B.mrowind.end() );
  stoprow = max_element( B.mrowind.begin(), B.mrowind.end() );

  // To create the A.mcolptr vector we loop through the existing row
  // numbers of B and count how many there are in B.mrowind.
  // First we make sure it is initialized to all zeros.
  A.mcolptr.assign( A.mcr+1, 0 );
  for (Index i=*startrow; i<=*stoprow; i++)
    {
      Index n = count( B.mrowind.begin(), B.mrowind.end(), i);

      // The column pointer for the column above the current one in A
      // are then set to the current one plus this amount.
      A.mcolptr[i+1] = A.mcolptr[i] + n;
    }

  // Next we loop through the columns of A and search for the corresponding
  // row index within each column of B (i.e. two loops). The elements are
  // then stored in A. We know that for every column in B we will have
  // the corresponding rowindex in A.
  vector<Numeric>::iterator dataA_it = A.mdata.begin();
  vector<Index>::iterator rowindA_it = A.mrowind.begin();

  // Looping over columns in A to keep track of what row number we are
  // looking for.
  for (size_t c=0; c<A.mcolptr.size(); ++c)
    {
      // Looping over columns in B to get the elements in the right order,
      // this will turn into row index in A.
      for (size_t i=0; i<B.mcolptr.size()-1; i++)
        {
          // Since only one element can occupy one row in a column, we
          // only need to call find() once per column in B.
          vector<Index>::const_iterator elem_it =
            find(B.mrowind.begin()+*(B.mcolptr.begin()+i),
                 B.mrowind.begin()+*(B.mcolptr.begin()+i+1), Index (c));

          // If we found the element, store it in A.mrowind and A.mdata
          if (elem_it != B.mrowind.begin()+*(B.mcolptr.begin()+i+1))
            {
              *rowindA_it = i;
              rowindA_it++;

              // To get the corresponding element in B.mdata we subtract
              // the initial value of the row index iterator from elem_it
              // and add the difference to a mdata iterator.
              Index diff = elem_it - B.mrowind.begin();
              vector<Numeric>::const_iterator elemdata_it=B.mdata.begin() + diff;
              *dataA_it = *elemdata_it;
              dataA_it++;
            }
        }
    }
}


//! Sparse - Sparse multiplication.
/*!
  Calculates A = B*C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix.
  \param B First product matrix.
  \param C Second product matrix.

  \author Mattias Ekstroem
  \date   2003-08-06
*/
void mult( Sparse& A,
           const Sparse& B,
           const Sparse& C )
{
  // Check dimensions and make sure that A is empty
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );
  A.mcolptr.assign( A.mcr+1, 0);
  A.mrowind.clear();
  A.mdata.clear();

  // Transpose B to simplify multiplication algorithm, after transposing we
  // can extract columns form the two matrices and multiply them, (which is
  // easier than extacting rows.)
  Sparse Bt(B.ncols(), B.nrows());
  //  cout << "B.nrows, B.ncols = " << B.nrows() << ", " << B.ncols() << "\n";
  //  cout << "B = \n" << B << "\n";
  transpose(Bt,B);

  // By looping over columns in C and multiply them with every column in Bt
  // (instead of the conventional loooping over Bt), we get the output
  // elements in the right order for storing them.
  for (size_t c=0; c<C.mcolptr.size()-1; ++c)
    {
      // Get row indices of this column
      //Index beginC = (*C.mcolptr)[c];
      //Index endC = (*C.mcolptr)[c+1];

      for (size_t b=0; b<Bt.mcolptr.size()-1; ++b)
        {
          // Get the intersection between the elements in the two columns and
          // and store them in a temporary vector
          std::set<Index> colintersec;
          set_intersection(C.mrowind.begin()+*(C.mcolptr.begin()+c),
            C.mrowind.begin()+*(C.mcolptr.begin()+c+1),
            Bt.mrowind.begin()+*(Bt.mcolptr.begin()+b),
            Bt.mrowind.begin()+*(Bt.mcolptr.begin()+b+1),
            inserter(colintersec, colintersec.begin()));

          // If we got an intersection, loop through it and multiply the
          // element pairs from C and Bt and store result in A
          if (!colintersec.empty())
            {
              Numeric tempA = 0.0;
              for (std::set<Index>::iterator i=colintersec.begin();
                i!=colintersec.end(); ++i)
                {
                  // To get iterators to the data elements in Bt and C, we
                  // subtract the iterator that points to the start of the
                  // mrowind vector from the iterators that points to the
                  // values in colintersec. We then get the number of steps 
                  // from the start of the both vectors mrowind and mdata to
                  // the intersecting values and can add them to the
                  // iterator that points to the beginning of mdata.
                  vector<Index>::const_iterator rowindBt_it =
                    find(Bt.mrowind.begin()+*(Bt.mcolptr.begin()+b),
                      Bt.mrowind.begin()+*(Bt.mcolptr.begin()+b+1), *i);
                  vector<Index>::const_iterator rowindC_it =
                    find(C.mrowind.begin()+*(C.mcolptr.begin()+c),
                      C.mrowind.begin()+*(C.mcolptr.begin()+c+1), *i);

                  vector<Numeric>::const_iterator dataBt_it =
                    Bt.mdata.begin()+(rowindBt_it-Bt.mrowind.begin());
                  vector<Numeric>::const_iterator dataC_it =
                    C.mdata.begin()+(rowindC_it-C.mrowind.begin());

                  tempA += *dataBt_it * *dataC_it;
                }
              A.rw(b,c) = tempA;
            }
        }
    }
}


//! Sparse - Sparse addition.
/*!
  Calculates A = B+C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of B, and C must match.  A will be resized.

  \param A Output: Result matrix.
  \param B First summand matrix.
  \param C Second summand matrix.

  \author Oliver Lemke
  \date   2009-09-03
*/
void add( Sparse& A,
          const Sparse& B,
          const Sparse& C )
{
  // Check dimensions
  assert( B.ncols() == C.ncols() );
  assert( B.nrows() == C.nrows() );

  // We copy the smaller matrix of B or C to A. This way we can
  // loop over the matrix with the fewer number of elements to perform
  // the actual addition later. 
  const Sparse* D;
  if (B.data().size() < C.data().size())
    {
      A=C;
      D=&B;
    }
  else
    {
      A=B;
      D=&C;
    }

  for (size_t c = 0; c < D->mcolptr.size()-1; ++c)
    {
      // Loop through the elements in this column:
      for (Index i = D->mcolptr[c]; i < D->mcolptr[c+1]; ++i)
        {
          A.rw(D->mrowind[i], c) += D->mdata[i];
        }
    }
}


//! Sparse - Sparse subtraction.
/*!
  Calculates A = B-C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of B, and C must match.  A will be resized.

  \param A Output: Result matrix.
  \param B First subtrahend matrix.
  \param C Second subtrahend matrix.

  \author Oliver Lemke
  \date   2009-09-03
*/
void sub( Sparse& A,
          const Sparse& B,
          const Sparse& C )
{
  // Check dimensions
  assert( B.ncols() == C.ncols() );
  assert( B.nrows() == C.nrows() );

  A=B;

  for (size_t c = 0; c < C.mcolptr.size()-1; ++c)
    {
      // Loop through the elements in this column:
      for (Index i = C.mcolptr[c]; i < C.mcolptr[c+1]; ++i)
        {
          A.rw(C.mrowind[i], c) -= C.mdata[i];
        }
    }
}


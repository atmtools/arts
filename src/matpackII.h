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

#ifndef matpackII_h
#define matpackII_h

#include <vector>
#include "matpackI.h"
#include "mystring.h"

// Declare existance of some classes:
class bifstream;
class bofstream;

/*
  Notes:

  There are two different ways to index: 
  S.rw(3,4) = 1;                // Read and write
  cout << S.ro(3,4);            // Read only

  This distinction is necessary, because rw() creates elements if they
  don't already exist.

 */


/** The implementation of a sparse matrix. 

    Class Sparse is derived from this one. It also allocates storage. 

    The chosen storage format is the `compressed column' format. This
    is the same format used by Matlab. See Matlab User Guide for
    a description.
*/
class SparseView {
public:
  // Member functions:
  Index nrows() const;
  Index ncols() const;
  Index nnz()   const;

  // Index Operators:
  Numeric& rw(Index r, Index c);
  Numeric  ro(Index r, Index c) const;

  Numeric  operator() (Index r, Index c) const;

  // Friends:
  friend std::ostream& operator<<(std::ostream& os, const SparseView& v);
  friend void mult (VectorView y, const SparseView& M, const ConstVectorView& x );
  friend void mult (MatrixView A, const SparseView B, const MatrixView C );
  friend void mult (SparseView A, const SparseView B, const SparseView C );
  friend void transpose (SparseView A, const SparseView B );
  friend void transpose2 (SparseView A, const SparseView B );
  // IO functions must be friends:
  friend void xml_write_to_stream (ostream& os_xml, const Sparse& sparse, 
                                   bofstream *pbofs, const String &name);

protected:
  // Constructors:
  SparseView();
  SparseView(std::vector<Numeric> *data,
                  std::vector<Index>   *rowind,
                  std::vector<Index>   *colptr,
                  const Range&         rr,
                  const Range&         cr                 );

  // Data members:
  /** The actual data values. */
  std::vector<Numeric> *mdata;
  /** Row indices. */
  std::vector<Index> *mrowind;
  /** Pointers to first data element for each column. */
  std::vector<Index> *mcolptr;
  /** Range of rows. */
  Range mrr;
  /** Range of columns. */
  Range mcr;
};

/** The Sparse class. This is a SparseView that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from SparseView. Just the constructors and
    the destructor have to be defined separately here, since these are not
    inherited. */
class Sparse : public SparseView {
public:
  // Constructors:
  Sparse();
  Sparse(Index r, Index c);
  Sparse(const Sparse& m);

  // Resize function:
  void resize(Index r, Index c);

  // Destructor:
  ~Sparse();

  Numeric operator() (Index r, Index c) const
    { return SparseView::operator() (r, c); }
};


// Functions for general matrix operations
void mult( VectorView y,
           const SparseView& M,
           const ConstVectorView& x );

void mult( MatrixView A,
           const SparseView B,
           const MatrixView C );

void mult( SparseView A,
           const SparseView B,
           const SparseView C );

void transpose( SparseView A,
                const SparseView B );

#endif

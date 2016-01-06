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

/*!
  \file   test_utils.h
  \author Simon Pfreundschuh <simonpf@chlamers.se>
  \date   Sun May  3 20:47:42 2015

  \brief Utility functions for testing.
*/

#ifndef test_utils_h
#define test_utils_h

#include "matpackI.h"
#include <stdlib.h>
#include <time.h>

/** Random number class.

    Uses rand() to generate a pseudo-random integer and converts it to
    rand_type and maps it to the range [lo, hi]. The current calendar
    time at construction is used to seed the generator.

*/
template <class rand_type> class Rand
{

public:
Rand( rand_type lo,
      rand_type hi ) : low( lo ), range( hi - lo )
    {
        srand( rand() );
    }

    rand_type operator ()() const
    {
        rand_type r = (rand_type) ( ((Numeric) rand()) /
                                    ((Numeric) RAND_MAX) * (Numeric) range );
        return  low + r;
    }

/** Random Index class.

    Template specialization for values of type Index to avoid rounding
    problems.
*/
private:
    rand_type low, range;

};

template <> class Rand<Index>
{

public:
Rand( Index lo,
      Index hi ) : low( lo ), range( hi - lo )
    {
        // Avoid negative ranges.
        if (hi <= lo)
            range = 0;
        srand( rand() );
    }

    Index operator ()() const
    {
        return  low + rand() % (range + 1);
    }

private:
    Index low, range;

};

// Add noise to vector.
void add_noise( VectorView v,
                Numeric range );

// Fill matrix with random values.
void random_fill_matrix( MatrixView A,
                         Numeric range,
                         bool positive );

// Fill sparse matrix with random values.
void random_fill_matrix( Sparse &A,
                         Numeric range,
                         bool positive );

// Fill a dense and a sparse matrix with the identical, random values.
void random_fill_matrix( Matrix& A,
                         Sparse& B,
                         Numeric range,
                         bool positive );

// Fill matrix with random values symmetrically.
void random_fill_matrix_symmetric( MatrixView A,
                                   Numeric range,
                                   bool positive );

// Generate random, positive semi-definite matrix.
void random_fill_matrix_pos_def( MatrixView A,
                                  Numeric range,
                                  bool positive );

// Generate random, positive semi-definite matrix.
void random_fill_matrix_pos_semi_def( MatrixView A,
                                      Numeric range,
                                      bool positive );

// Fill vector with random values.
void random_fill_vector( VectorView A,
                         Numeric range,
                         bool positive );


// Pick random submatrix.
MatrixView random_submatrix( MatrixView A,
                             Index m,
                             Index n );

// Generate random range in the range [0, n - 1]
Range random_range( Index n );

// Maximum element-wise error of two matrices.
Numeric get_maximum_error( ConstMatrixView A1,
                           ConstMatrixView A2,
                           bool relative );

// Maximum element-wise error of two matrices.
Numeric get_maximum_error( ConstVectorView v1,
                           ConstVectorView v2,
                           bool relative );

#endif // test_utils_h

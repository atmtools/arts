/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
                      Patrick Eriksson <patrick@rss.chalmers.se>

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   vecmat.h

*/



#ifndef vecmat_h
#define vecmat_h

////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

//--------------------< Global ARTS header >--------------------
// Defines for example the type Numeric.
#include "arts.h"

//--------------------< STL headers: >--------------------
#include <vector>		// Vektor class from the STL

//--------------------< TNT headers: >--------------------

// Set the type of TNT indices to size_t, for compatibility with STL
// vectors and strings.
#define TNT_SUBSCRIPT_TYPE size_t

#include "tnt.h"
#include "vec.h"
#include "fmat.h"
#include "fcscmat.h"
// TNT stopwatch: (seems not to work anymore)
//#include "stpwatch.h"
// TNT Vector addapter for Array:
#include "vecadaptor.h"


using namespace TNT;

// Other important TNT types:
// Index1D: Can be used to select sub-ranges in matrices. 

// If NDEBUG is defined then we do not want to have any
// TNT bounds checks:
#ifdef NDEBUG
#define TNT_NO_BOUNDS_CHECK
#endif  



////////////////////////////////////////////////////////////////////////////
//   Some definitions
////////////////////////////////////////////////////////////////////////////

/** For numeric matrices. 
    Awailable matrix/vector classes are:
    \begin{itemize}
    \item MATRIX
    \item VECTOR
    \item ARRAY
    \end{itemize} 

    VECTORS can be efficiently multiplied by a MATRIX or by another
    VECTOR. They are made up of elements of type Numeric. ARRAY can be
    used in exactly the same way as VECTOR (including optional 1-based
    indexing and bounds checking), but is just meant to store
    things. Use VECTOR for numerics, ARRAY for all other vectors.
    @see VECTOR ARRAY */
typedef TNT::Fortran_Matrix<Numeric> MATRIX;

/** Sparse Matrix.
    @author Stefan Buheler
    @date   2000-08-16. */
//typedef Fortran_Sparse_Col_Matrix<Numeric> SPARSEMATRIX;


/** Regions for matrices.
    @author Stefan Buehler 11.06.2000. */
typedef TNT::Region2D<MATRIX> REGION2D;

/** Regions for const matrices.
    @author Stefan Buehler 11.06.2000. */
typedef TNT::const_Region2D<MATRIX> const_REGION2D;

/** For numeric vectors.
    @see MATRIX ARRAY */
typedef TNT::Vector<Numeric> VECTOR;

/** Regions for vectors.
    @author Stefan Buehler 11.06.2000. */
typedef TNT::Region2D<VECTOR> REGION1D;

/** Regions for const vectors.
    @author Stefan Buehler 11.06.2000. */
typedef TNT::const_Region2D<VECTOR> const_REGION1D;

/** Subscript type for MATRIX, VECTOR, and ARRAY.
    @author Stefan Buehler
    @date 2000-08-16 */
typedef TNT::Subscript SUBSCRIPT;

/** For arrays. This can be used in the same way as VECTOR, but can
    store arbitrary elements. Furthermore, this has the
    member function get_vector, which returns the underlying
    std::vector. Thus, all STL functions and algorithms can be used on
    this. WARNING! You should know what you are doing if you want to
    use STL functions directly. Some of them will not work or cause
    very strange results. (Basically all functions that re-allocate
    the memory so that pointers become invalid.)
    
    Because constructors are not inherited, I have to re-define all
    costructors.  

    @see MATRIX ARRAY */
template<class EE>
class ARRAY : public TNT::Vector_Adaptor< std::vector<EE> >
{
public:

  typedef std::vector<EE> mybase;

  ARRAY() : TNT::Vector_Adaptor<mybase>() {};
  ARRAY(const ARRAY<EE> &A) :
    TNT::Vector_Adaptor<mybase>(A) {};
  
  ARRAY(SUBSCRIPT N, /*const*/ char *s) :
    TNT::Vector_Adaptor<mybase>(N, s) {};
  
  ARRAY(SUBSCRIPT N, const EE& value = EE()) :
    TNT::Vector_Adaptor<mybase>(N, value) {};

  ARRAY(SUBSCRIPT N, const EE* values) :
    TNT::Vector_Adaptor<mybase>(N, values) {};

  ARRAY(const EE & A) :
    TNT::Vector_Adaptor<mybase>(A) {};
};



////////////////////////////////////////////////////////////////////////////
//   TNT patches
//
//   This part includes patches for TNT to fulfill the functionality
//   stated above. 
//
//   2000-04-11 Patrick Eriksson. Adaption of earlier version.
//   2000-08-16 Stefan Buehler: 
//      Converted template functions to direct functions of type MATRIX
//      and VECTOR.   
//
////////////////////////////////////////////////////////////////////////////

//// MATRIX - MATRIX ///////////////////////////////////////////////////////

MATRIX emult(const MATRIX &A, const MATRIX &B);

MATRIX ediv(const MATRIX &A, const MATRIX &B);



//// VECTOR - VECTOR ///////////////////////////////////////////////////////

VECTOR emult( const VECTOR &A, const VECTOR &B);

VECTOR ediv( const VECTOR &A, const VECTOR &B);



//// MATRIX - SCALAR ///////////////////////////////////////////////////////

MATRIX operator+(const MATRIX &A, const Numeric scalar);

MATRIX operator+(const Numeric scalar, const MATRIX &A);

MATRIX operator-(const MATRIX &A, const Numeric scalar);

MATRIX operator-(const Numeric scalar, const MATRIX &A);

MATRIX operator*(const MATRIX &A, const Numeric scalar);

MATRIX operator*(const Numeric scalar, const MATRIX &A);

MATRIX operator/(const MATRIX &A, const Numeric scalar);

MATRIX operator/(const Numeric scalar, const MATRIX &A);
        


//// VECTOR - SCALAR ///////////////////////////////////////////////////////

VECTOR operator+(const VECTOR &A, const Numeric scalar);

VECTOR operator+(const Numeric scalar, const VECTOR &A);

VECTOR operator-(const VECTOR &A, const Numeric scalar);

VECTOR operator-(const Numeric scalar, const VECTOR &A);

VECTOR operator*(const VECTOR &A, const Numeric scalar);

VECTOR operator*(const Numeric scalar, const VECTOR &A);

VECTOR operator/(const VECTOR &A, const Numeric scalar);

VECTOR operator/(const Numeric scalar, const VECTOR &A);



////////////////////////////////////////////////////////////////////////////
//   Definition of arrays of matrices and vectors
////////////////////////////////////////////////////////////////////////////

/** An array of matrices. */
typedef ARRAY<MATRIX> ARRAYofMATRIX;

/** An array of vectors. */
typedef ARRAY<VECTOR> ARRAYofVECTOR;

/** An array of integers. */
typedef ARRAY<size_t> ARRAYofsizet;

/** An array of strings. */
typedef ARRAY<string> ARRAYofstring;




#endif // vecmat_h

#ifndef vecmat_h
#define vecmat_h

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
#include "cmat.h"
// TNT stopwatch: (seems not to work anymore)
//#include "stpwatch.h"
// TNT Vector addapter for Array:
#include "vecadaptor.h"


// using namespace TNT;

// Other important TNT types:
// Index1D: Can be used to select sub-ranges in matrices. 

// If NDEBUG is defined then we do not want to have any
// TNT bounds checks:
#ifdef NDEBUG
#define TNT_NO_BOUNDS_CHECK
#endif  

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
typedef TNT::Matrix<Numeric> MATRIX;

/** For numeric vectors.
    @see MATRIX ARRAY */
typedef TNT::Vector<Numeric> VECTOR;

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
  
  ARRAY(TNT::Subscript N, /*const*/ char *s) :
    TNT::Vector_Adaptor<mybase>(N, s) {};
  
  ARRAY(TNT::Subscript N, const EE& value = EE()) :
    TNT::Vector_Adaptor<mybase>(N, value) {};

  ARRAY(TNT::Subscript N, const EE* values) :
    TNT::Vector_Adaptor<mybase>(N, values) {};

  ARRAY(const EE & A) :
    TNT::Vector_Adaptor<mybase>(A) {};
};


#endif // vecmat_h

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


using namespace TNT;

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


/*------------------< Allowed MATRIX and VECTOR operations >------------------

  This part defines the MATRIX and VECTOR functionality that will be 
  maintained in ARTS. Other features can work, but avoid these as they
  can be removed in future versions.

  The following definitions are used below
    MATRIX   A, B, C;
    VECTOR   x, y, z;
    Numeric  s;
    size_t   i, j;
    Index1D  I, J, K, L;


  Initialization and resizing.

    MATRIX A(i,j)
             Creates a matrix with i rows and j columns. The values of A
             are undefined.    

    MATRIX A(i,j,s)
             Creates a matrix with i rows and j columns, and sets all the
             elements of A to s.

    MATRIX A(2, 4, " 1  2  0  4 "
                   " 2  0  9  7 ");
             This example shows how to create a matrix with 2 rows and 4
             columns with specified element values.

    A.newsize(i,j)
             Resizes the matrix to have i rows and j columns. The values of
             A are undefined.    

    B = A
             Makes B a copy of A, including resizing of B.

    A = s
             Sets all elements of A to s.

    VECTOR x(i)
             Creates a vector of length i. The values of x are undefined.    

    VECTOR x(i,s)
             Creates a matrix of length i, and sets all the elements of 
             x to s.

    VECTOR x(4, " 1  2  0  4 " );
             This example shows how to create a vector of length 4 with
             specified element values.

    x.newsize(i)
             Resizes the vector to have length i. The values of x are 
             undefined.    

    y = x
             Makes y a copy of x, including resizing of y.

    x = s
             Sets all elements of x to s.


  Dimensions

    i = A.dim(1), j = A.dim(2)
             The size of a matrix. The dimension order is row - columns, 
             i.e. i is the number of rows and j is the number of columns. 

    i = x.dim()
             The length of a vector.


  Basic math:

    C = A*B  
             Matrix multiplication. The number of columns of A must be equal
             to the number of rows of B.

    C = A+B, A-C  
             Summation or difference of matrices. The both matrices must
             have the same size.

    C = emult(A,B)
             Element-by-element product of matrices, i.e. C(i,j) = 
             A(i,j)*B(i,j). A and B must have the same size.

    C = ediv(A,B)
             Element-by-element division of matrices, i.e. C(i,j) = 
             A(i,j)/B(i,j). A and B must have the same size.

    C = transpose(A)
             Transpose of a matrix.

    y = A*x  
             Matrix vector multiplication. The number of columns of A must be 
             equal to the length of x. The reversed order, x*A, is not allowed.

    z = x+y, x-y  
             Summation or difference of vectors. The two vectors must
             have the same length.

    z = emult(x,y)
             Element-by-element product of vectors, i.e. z(i) = x(i)*y(i).
             The vectors must have the same length. 

    z = ediv(x,y)
             Element-by-element division of vectors, i.e. z(i) = x(i)/y(i) 
             The vectors must have the same length. 

    s = dot_prod(x,y)   
             Scalar product of two vectors. The vectors must have the same 
             length. The result is s=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)+... 
    
    C = s+A, A+s, s-A, A-s, s*A, A*s, s/A or A/s
             Summation, difference, multiplication or division with a scalar
             of each element of a matrix. 

    y = s+x, x+s, s-x, x-s, s*x, x*s, s/x or x/s
             Summation, difference, multiplication or division with a scalar
             of each element of a vector. 


  Conversion between vectors and matrices:

    x = to_vector(A), to_vector(x,A)
             Converts a matrix to a vector. The matrix can either be a 
             column (n x 1) or row (1 x n) vector.

    A = to_matrix(x), A = to_matrix(x)
             Converts a vector of length n to a matrix of size (n x 1),
             i.e. the vector is interpreted as a column.


  Index and regions
    Vector and matrix indecies are put between parenthesis. The first 
    element has index 1 (not 0).

    s = A(i,j)  
             Gives element j of row i in A.    

    C(I,J) = A(K,L)
             Picks out a matrix region (with natural limitations regarding
             the size of A and C). The matrix on the left hand side must be
             indexed ( i.e. C=A(K,L) is not valid even if sizes match).

    A(I,J+2), A(I-2,J)
             An integer value can be added to or removed from the index 
             ranges ((with natural limitations regarding the size of A).
             
    A(I,J) = s
             A matrix region can be set to a value.

    s = x(i)     
             Gives element i of the vector x.


  11.04.00 Written by Patrick Eriksson.
---------------------------------------------------------------------------*/



/*------------------------------< TNT PATCHES >---------------------------*/

// This part includes patches for TNT to fulfill the functionality
// stated above. 
//
// 11.04.00 Patrick Eriksson. Adaption of earlier version.

// BASIC MATH
//--------------------------------------------------------------------------

// MATRIX - MATRIX ---------------------------------------------------------

// emult
//
  template <class T>
  Matrix<T> emult(const Matrix<T> &A, const Matrix<T> &B)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    assert(M==B.dim(1));
    assert(N==B.dim(2));

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] * B[i][j];

    return tmp;
  }        

// ediv
//
  template <class T>
  Matrix<T> ediv(const Matrix<T> &A, const Matrix<T> &B)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    assert(M==B.dim(1));
    assert(N==B.dim(2));

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] / B[i][j];

    return tmp;
  }        


// VECTOR - VECTOR ---------------------------------------------------------

// emult
//
  template <class T>
  Vector<T> emult( const Vector<T> &A, const Vector<T> &B)
  {
    Subscript N = A.dim();

    assert(N==B.dim());

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
            tmp[i] = A[i] * B[i];

    return tmp;
  }       

// ediv
//
  template <class T>
  Vector<T> ediv( const Vector<T> &A, const Vector<T> &B)
  {
    Subscript N = A.dim();

    assert(N==B.dim());

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
            tmp[i] = A[i] / B[i];

    return tmp;
  }       


// MATRIX - SCALAR ---------------------------------------------------------

// + 
//
  template <class T>
  Matrix<T> operator+(const Matrix<T> &A, const T& scalar)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] + scalar;

    return tmp;
  }
  template <class T>
  Matrix<T> operator+(const T& scalar, const Matrix<T> &A)
  {
    return ( A + scalar );
  }                 

// - 
//
  template <class T>
  Matrix<T> operator-(const Matrix<T> &A, const T& scalar)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] - scalar;

    return tmp;
  }
  template <class T>
  Matrix<T> operator-(const T& scalar, const Matrix<T> &A)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] =  scalar - A[i][j];

    return tmp;
  }


// *
//
  template <class T>
  Matrix<T> operator*(const Matrix<T> &A, const T& scalar)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] * scalar;

    return tmp;
  }
  template <class T>
  Matrix<T> operator*(const T& scalar, const Matrix<T> &A)
  {
    return ( A * scalar );
  }                 

//
// /
  template <class T>
  Matrix<T> operator/(const Matrix<T> &A, const T& scalar)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] / scalar;

    return tmp;
  }                       
  template <class T>
  Matrix<T> operator/(const T& scalar, const Matrix<T> &A)
  {
    Subscript M = A.dim(1);
    Subscript N = A.dim(2);

    Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = scalar / A[i][j];

    return tmp; 
  }        


// VECTOR - SCALAR ---------------------------------------------------------

// + 
//
  template <class T>
  Vector<T> operator+(const Vector<T> &A, const T& scalar)
  {
    Subscript N = A.dim();

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
      tmp[i] = A[i] + scalar;

    return tmp;
  }
  template <class T>
  Vector<T> operator+(const T& scalar, const Vector<T> &A)
  {
    return ( A + scalar );
  }

// -
//
  template <class T>
  Vector<T> operator-(const Vector<T> &A, const T& scalar)
  {
    Subscript N = A.dim();

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
      tmp[i] = A[i] - scalar;

    return tmp;
  }
  template <class T>
  Vector<T> operator-(const T& scalar, const Vector<T> &A)
  {
    Subscript N = A.dim();

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
      tmp[i] = scalar - A[i];

    return tmp;
  }

// *
//
  template <class T>
  Vector<T> operator*(const Vector<T> &A, const T& scalar)
  {
    Subscript N = A.dim();

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
      tmp[i] = A[i] * scalar;

    return tmp;
  }
  template <class T>
  Vector<T> operator*(const T& scalar, const Vector<T> &A)
  {
    return ( A * scalar );
  }

// /
//
  template <class T>
  Vector<T> operator/(const Vector<T> &A, const T& scalar)
  {
    Subscript N = A.dim();

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
      tmp[i] = A[i] / scalar;

    return tmp;
  }
  template <class T>
  Vector<T> operator/(const T& scalar, const Vector<T> &A)
  {
    Subscript N = A.dim();

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
      tmp[i] = scalar / A[i];

    return tmp;
  }



// MATRIX AND VECTOR CONVERSIONS
//--------------------------------------------------------------------------
//
// The functions to_matrix and to_vector are placed in math_funcs.cc.


/*------------------< End MATRIX and VECTOR operations >----------------*/


/*---------------< Define Arrays of Matrices and Vectors >---------------*/

/** An array of matrices. */
typedef ARRAY<MATRIX> ARRAYofMATRIX;

/** An array of vectors. */
typedef ARRAY<VECTOR> ARRAYofVECTOR;

// You can write variables of these types to stdout, simply by 
// cout << x;



#endif // vecmat_h

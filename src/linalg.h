#ifndef linalg_h
#define linalg_h

#include "mtl/mtl.h"
#include "mtl/utils.h"
#include "mtl/matrix.h"
#include "mtl/linalg_vec.h"

//using namespace mtl;
//using namespace std;

// typedef mtl::matrix< Numeric, 
//   rectangle<>, 
//   dense<>, 
//   row_major>::type MATRIX; 

typedef mtl::dense1D<Numeric> VEC;

typedef mtl::matrix<Numeric>::type MAT;

typedef mtl::matrix<Numeric,mtl::rectangle<>,mtl::compressed<> >::type SPARSE_MAT;

#define ARRAY mtl::dense1D       


/** Output function for vectors (for example VEC and ARRAY).

    This function is taken directly from MTL, with the only difference
    that it can print to any output stream, not just to stdout.  Only
    minor adaptations to ARTS coding standards have been applied ("\n"
    instead of endl).

    There is a comma after the last element of the vector, because
    this makes it simpler. Could be changed for nicer optics later. 

    \param os The output stream.
    \param x  The vector to print.

    \author Stefan Buehler
    \date   2000-12-06
 **/
template <class Vector>
inline void
print_vector(ostream& os, Vector x)
{
  using namespace mtl;
  typename Vector::iterator i = x.begin();
  std::cout << "[";
  while (not_at(i, x.end())) {
    std::cout << *i << ",";
    ++i;
  }
  std::cout << "]" << std::endl;
}

/** Input function for vectors (for example VEC and ARRAY).

    This functions reads an MTL vector from an input stream. The format is
    very simple, first there should be an integer indicating the
    length of the vector, then the elements. Thus, it will work for
    arbitrary element types, providing their input operator is
    overloaded correctly. 

    \param is The input stream.
    \retval x  The vector to read.

    \author Stefan Buehler
    \date   2000-12-08
 **/
template <class Vector>
inline void
read_vector_from_stream(Vector x, istream& is)
{
  INDEX n;
  is >> n;
  x.resize(n);
  for(INDEX i=0; i<n; ++i)
    {
      is >> x[i];
    }
}


/** Erase an element of a vector (for example VEC and ARRAY). This
    functionality is missing in MTL vectors, but we need in some
    places. 


    \retval x  The vector to treat.
    \param k   Index of the element to erase.

    \author Stefan Buehler
    \date   2000-12-08
 **/
template <class Vector>
inline void
erase_vector_element(Vector& x, INDEX k)
{
  Vector y(x.size()-1);
  copy(x(0,k),y(0,k));
  copy(x(k+1,x.size()),y(k,x.size()-1));
  x = y;
}


/** Output function for matrices.

    This function is taken directly from MTL, with the only difference
    that it can print to any output stream, not just to stdout.  Only
    minor adaptations to ARTS coding standards have been applied ("\n"
    instead of endl).

    \param os The output stream.
    \param A  The matrix to print.

    \author Stefan Buehler
    \date   2000-12-06
 **/
template <class Matrix>
inline void
print_all_matrix(ostream& os, const Matrix& A)
{
  using namespace mtl;
  typedef typename matrix_traits<Matrix>::size_type Int;
  Int i,j;
  std::cout << A.nrows() << "x" << A.ncols() << std::endl;
  std::cout << "[" << std::endl;
  for (i=0; i < A.nrows(); ++i) {
    std::cout << "[";
    for (j=0; j < A.ncols(); ++j) {
      std::cout << A(i,j);
      if (j < A.ncols() - 1)
	std::cout << ",";
    }
    std::cout << "]";
    if (i < A.nrows() - 1)
      std::cout << "," << std::endl;
    else
      std::cout << std::endl;
  }
  std::cout << "]" << std::endl;
}




#endif  // linalg_h

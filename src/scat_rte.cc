//! Calculate vector radiative transfer with fixed scattering integral
/*! 
  This function computes the radiative transfer for a thin layer.
  It is a general function which works for both, the vector and the scalar 
  RTE. But for the scalar equation it is more efficient to use the method
  *stokes_vecScalar*.
  All coefficients and the scattered field vector are assumed to be constant
  inside the grid cell/layer.
  Then an analytic solution can be found (see AUG for details). 
  
  \param stokes_vec Output and Input: Stokes Vector after traversing a grid
                 cell/layer. 
  \param ext_mat Input: Extinction coefficient matrix.
  \param abs_vec Input: Absorption coefficient vector.
  \param sca_vec Input: Scattered field vector.
  \param l_step  Input: Pathlength through a grid cell/ layer.
  \param a_planck_value  Input: Planck function.
  \param stokes_dim Input: Stokes dimension. 

  \author Claudia Emde
  \date 2002-06-08
*/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "logic.h"
#include "lin_alg.h"


void
stokes_vecGeneral(//WS Output and Input:
               Vector& stokes_vec,
               //WS Input:
               const Matrix& ext_mat,
               const Vector& abs_vec,
               const Vector& sca_vec,
               const Numeric& l_step,
               const Numeric& a_planck_value,
               const Index& stokes_dim)
{ 

  // Check if stokes_vec has stokes_dim components
  assert(is_size(stokes_vec, stokes_dim));
  // check, if ext_mat is quadratic
  assert(is_size(ext_mat, stokes_dim, stokes_dim)); 
  // check if the dimensions agree
  assert(is_size(abs_vec, stokes_dim));
  assert(is_size(sca_vec, stokes_dim));

  //Initialize internal variables:

  // Matrix LU used for LU decompostion and as dummy variable:
  Matrix LU(stokes_dim, stokes_dim); 
  ArrayOfIndex indx(stokes_dim); // index for pivoting information 
  Vector b(stokes_dim); // dummy variable 
  Vector x(stokes_dim); // solution vector for K^(-1)*b
  Matrix I(stokes_dim, stokes_dim);

  Vector B_abs_vec(stokes_dim);
  B_abs_vec = abs_vec;
  B_abs_vec *= a_planck_value; 
  
  for (Index i=0; i<stokes_dim; i++) 
    b[i] = B_abs_vec[i] + sca_vec[i];  // b = abs_vec * B + sca_vec

  // solve K^(-1)*b = x
  ludcmp(LU, indx, ext_mat);
  lubacksub(x, LU, b, indx);

  Matrix ext_mat_ds(stokes_dim, stokes_dim);
  ext_mat_ds = ext_mat;
  ext_mat_ds *= -l_step; // ext_mat = -ext_mat*ds
  
  Index q = 10;  // index for the precision of the matrix exponential function
  Matrix exp_ext_mat(stokes_dim, stokes_dim);
  matrix_exp(exp_ext_mat, ext_mat_ds, q);

  Vector term1(stokes_dim);
  Vector term2(stokes_dim);
  
  id_mat(I);
  for(Index i=0; i<stokes_dim; i++)
    {
      for(Index j=0; j<stokes_dim; j++)
	LU(i,j) = I(i,j) - exp_ext_mat(i,j); // take LU as dummy variable
    }
  mult(term2, LU, x); // term2: second term of the solution of the RTE with
                      //fixed scattered field

  
  mult(term1, exp_ext_mat, stokes_vec);// term1: first term of the solution of
                                    // the RTE with fixed scattered field
  
  for (Index i=0; i<stokes_dim; i++) 
    stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector

  
}



//! Calculate scalar radiative transfer with fixed scattering integral
/*! 
  This function computes the radiative transfer for a thin layer.
  All coefficients and the scattered field vector are assumed to be constant
  inside the grid cell/layer.
  Then an analytic solution can be found (see AUG for details). 
  
  \param stokes_vec Output and Input: Stokes Vector after traversing a grid
  cell/layer. 
  \param ext_mat Input: Extinction coefficient matrix.
  \param abs_vec Input: Absorption coefficient vector.
  \param sca_vec Input: Scattered field vector.
  \param l_step  Input: Pathlength through a grid cell/ layer.
  \param a_planck_value  Input: Planck function.
  \param stokes_dim Input: Stokes dimension.
  
  \author Claudia Emde
  \date 2002-06-08
*/
void
stokes_vecScalar(//WS Input and Output:
	      Vector& stokes_vec,
	      //WS Input: 
	      const Matrix& ext_mat,
	      const Vector& abs_vec,
	      const Vector& sca_vec,
	      const Numeric& l_step,
	      const Numeric& a_planck_value,
	      const Index& stokes_dim)
{ 
  //Check if we really consider the scalar case.
  if(stokes_dim != 1)
    throw runtime_error("You can use the method *stokes_vecScalar* only if"
                        "*stokes_dim* equals 1"); 

  // Check stokes_vec, one component vector:
  assert(is_size(stokes_vec, stokes_dim));
  // Check, if ext_mat is a scalar, i.e. 1x1 matrix
  assert(is_size(ext_mat, stokes_dim, stokes_dim)); 
  // Check, if  absorption coefficient and scattering coefficients   
  // coefficients are 1 component vectors.                       
  assert(is_size(abs_vec, stokes_dim)); 
  assert(is_size(sca_vec, stokes_dim));
  
  //Introduce scalar variables:
  
  
  //Extinction coefficient:
  Numeric ext_coeff = ext_mat(0,0);
  //Absorption coefficient:
  Numeric abs_coeff = abs_vec[0];
  //Scalar scattering integral:
  Numeric sca_int1D = sca_vec[0];

  //Intensity:
  Numeric stokes_vec1D = stokes_vec[0];
  
  //Do a radiative transfer step calculation:
  stokes_vec1D = stokes_vec1D * exp(-ext_coeff*l_step) + 
    (abs_coeff*a_planck_value+sca_int1D) /
    ext_coeff* (1-exp(-ext_coeff*l_step));
 
  //Put the first component back into *sto_vec*:
  stokes_vec[0] = stokes_vec1D;
  
}
		

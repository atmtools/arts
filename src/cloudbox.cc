/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
   USA. 
*/

/*!
  \file   cloudbox.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Thu May  23 10:59:55 2002
  
  \brief  Internal functions for scattering calculations.
*/

#include "cloudbox.h"

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <stdexcept>
#include <math.h>
#include "lin_alg.h"
#include "arts.h"
#include "auto_md.h"
#include "make_vector.h"
#include "array.h"
#include "logic.h"
#include "ppath.h"
#include "interpolation.h"
#include "physics_funcs.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/
//! Calculate vector radiative transfer with fixed scattering integral
/*! 
  This function computes the radiative transfer for a thin layer. All
  coefficients are assumed to be constant inside the grid cell/layer.
  For a fixed scattered field inside the layer an 
  analytic solution can be found (see AUG). 
  
  \param sto_vec Output and Input: Stokes Vector after traversing a grid
                 cell/layer. 
  \param ext_mat Input: Extinction coefficient matrix.
  \param abs_vec Input: Absorption coefficient vector.
  \param sca_vec Input: Scattered field vector.
  \param ds      Input: Pathlength through a grid cell/ layer.
  \param B       Input: Planck function.
  \param stokes_dim Input: Stokes dimension. 
*/
void
sto_vecCalc(VectorView sto_vec,
	    ConstMatrixView ext_mat,
	    ConstVectorView abs_vec,
	    ConstVectorView sca_vec,
	    const Numeric& ds,
	    const Numeric& B,
	    const Index& stokes_dim)
{ 
  // Stokes dimension
  assert(stokes_dim <= 4 && stokes_dim != 1);
 
  // check, if ext_mat is quadratic
  assert(is_size(ext_mat, stokes_dim, stokes_dim)); 
  // check if the dimensions agree
  assert(is_size(abs_vec, stokes_dim));
  assert(is_size(sca_vec, stokes_dim));


  //Initialize internal variables:

  // Matrix LU used for LU decompostion and as dummy variable:
  Matrix LU(stokes_dim,stokes_dim); 
  ArrayOfIndex indx(stokes_dim); // index for pivoting information 
  Vector b(stokes_dim); // dummy variable 
  Vector x(stokes_dim); // solution vector for K^(-1)*b
  Matrix I(stokes_dim,stokes_dim);

  Vector B_abs_vec(stokes_dim);
  B_abs_vec = abs_vec;
  B_abs_vec *= B; 
  
  for (Index i=0; i<stokes_dim; i++) 
    b[i] = abs_vec[i] + sca_vec[i];  // b = abs_vec * B + sca_vec

  // solve K^(-1)*b = x
  ludcmp(LU, indx, ext_mat);
  lubacksub(x, LU, b, indx);

  Matrix ext_mat_ds(stokes_dim,stokes_dim);
  ext_mat_ds = ext_mat;
  ext_mat_ds *= -ds; // ext_mat = -ext_mat*ds

  Index q = 10;  // index for the precision of the matrix exponential function
  Matrix exp_ext_mat(stokes_dim,stokes_dim);
  matrix_exp(exp_ext_mat, ext_mat, q);

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

  
  mult(term1, exp_ext_mat, sto_vec);// term1: first term of the solution of
                                    // the RTE with fixed scattered field
  
  for (Index i=0; i<stokes_dim; i++) 
    sto_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector

    
}



//! Calculate scalar radiative transfer with fixed scattering integral
/*! 
  This function computes the radiative transfer for a thin layer. All
  coefficients are assumed to be constant inside the grid cell/layer.
  For a fixed scattered field inside the layer an 
  analytic solution can be found (see AUG). 
  
  \param sto_vec Output and Input: Stokes Vector after traversing a grid
                 cell/layer. 
  \param ext_mat Input: Extinction coefficient matrix.
  \param abs_vec Input: Absorption coefficient vector.
  \param sca_vec Input: Scattered field vector.
  \param ds      Input: Pathlength through a grid cell/ layer.
  \param B       Input: Planck function.
  \param stokes_dim Input: Stokes dimension.
*/
void
sto_vec1DCalc(VectorView sto_vec,
	      ConstMatrixView ext_mat,
	      ConstVectorView abs_vec,
	      ConstVectorView sca_vec,
	      const Numeric& ds,
	      const Numeric& B,
	      const Index& stokes_dim)
{ 
  //Check if we really consider the scalar case.
  assert(stokes_dim == 1);

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
  Numeric sto_vec1D = sto_vec[0];
  
  //Do a radiative transfer step calculation:
  sto_vec1D = sto_vec1D * exp(-ext_coeff*ds) + (abs_coeff*B+sca_int1D) /
    ext_coeff* (1-exp(-ext_coeff*ds));
 
  //Put the first component back into *sto_vec*:
  sto_vec[0] = sto_vec1D;
}






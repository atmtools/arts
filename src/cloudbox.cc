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
*/
void rte_scat_vecCalc(VectorView sto_vec,
		      ConstMatrixView ext_mat,
		      ConstVectorView abs_vec,
		      ConstVectorView sca_vec,
		      const Numeric& ds,
		      const Numeric& B)
{ 
  // Stokes dimension
  const Index dim = ext_mat.nrows();
  assert(dim <= 4); 
  
  assert(is_size(ext_mat, dim, dim)); // check, if ext_mat is quadratic
  // check if the dimensions agree
  assert(is_size(abs_vec, dim));
  assert(is_size(sca_vec, dim));


  //Initialize internal variables
  Matrix LU(dim,dim); // used for LU decompostion and as dummy variable
  ArrayOfIndex indx(dim); // index for pivoting information 
  Vector b(dim); // dummy variable 
  Vector x(dim); // solution vector for K^(-1)*b
  Matrix I(dim,dim);

  Vector B_abs_vec(dim);
  B_abs_vec = abs_vec;
  B_abs_vec *= B; 
  
  for (Index i=0; i<dim; i++) 
    b[i] = abs_vec[i] + sca_vec[i];  // b = abs_vec * B + sca_vec

  // solve K^(-1)*b = x
  ludcmp(LU, indx, ext_mat);
  lubacksub(x, LU, b, indx);

  Matrix ext_mat_ds(dim,dim);
  ext_mat_ds = ext_mat;
  ext_mat_ds *= -ds; // ext_mat = -ext_mat*ds

  Index q = 10;  // index for the precision of the matrix exponential function
  Matrix exp_ext_mat(dim,dim);
  matrix_exp(exp_ext_mat, ext_mat, q);

  Vector term1(dim);
  Vector term2(dim);
  
  id_mat(I);
  for(Index i=0; i<dim; i++)
    {
      for(Index j=0; j<dim; j++)
	LU(i,j) = I(i,j) - exp_ext_mat(i,j); // take LU as dummy variable
    }
  mult(term2, LU, x); // term2: second term of the solution of the RTE with
                      //fixed scattered field

  
  mult(term1, exp_ext_mat, sto_vec);// term1: first term of the solution of
                                    // the RTE with fixed scattered field
  
  for (Index i=0; i<dim; i++) 
    sto_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector

    
}





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


//! 1D RT calculation inside the cloud box.
/*! 
  This function loops over all grid points and all directions and performs 
  the RT calculation with a fixed scattering integral for one frequency 
  of the frequency grid specified by *f_index*. 

  Note: The function uses the same input and output variables as the
  equivalent function for the 
  3D calculations. Here the dimensions which are not needed are empty.
  
  The agendas for computing ext_mat, pha_mat and abs_vec are needed as further
  input!
  
  \param i_field       Output: Updated intensity field. 
  \param i_field_old   Old intensity field.
  \param amp_mat       Amplitude matrix. 
  \param sca_field     Scattering integral at all grid points
                              inside the cloud box.
  \param cloudbox_limits Limits of the cloudbox.
  \param scat_za_grid  Zenith angle grid inside the cloud box.
  \param scat_aa_grid  Azimuthal angle grid inside the cloud box.
  \param p_grid        Pressure grid.
  \param lat_grid      Latitude grid.
  \param lon_grid      Longitude grid.
  \param t_field       Temperature field for alls grid points.
  \param z_field       Geometrical altitude field.
  \param z_ground      Ground altitude.
  \param r_geoid       Matrix containing geoids.
  \param f_grid        Frequency grid.
  \param f_index       Frequency index.
  \param blackbody_ground Flag to treat ground as blackbody.
  \param stokes_dim    The number of Stokes components to be calculated.
*/
void i_field_update1D(
		     Tensor6View i_field,
		     ConstTensor6View i_field_old,
		     ConstTensor7View amp_mat,
		     ConstTensor6View sca_field,
		     const ArrayOfIndex cloudbox_limits,
		     ConstVectorView scat_za_grid,
		     ConstVectorView scat_aa_grid,
		     ConstVectorView p_grid,
		     ConstVectorView lat_grid,
		     ConstVectorView lon_grid,
		     ConstTensor3View t_field,
		     ConstTensor3View z_field,
		     ConstMatrixView z_ground,
		     ConstMatrixView r_geoid,
		     ConstVectorView f_grid,
		     const Index f_index,
		     const Index blackbody_ground,
		     const Index stokes_dim
		     )
{

  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  

  //For the 1D Geometry these dimensions are empty.
  const Index scat_aa_index = 0;
  const Index lat_index = 0;
  const Index lon_index  = 0;


  //Loop over all directions, defined by scat_za_grid 
  for(Index scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      //Loop over all positions inside the cloud box defined by scat_p_grid
      for(Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	  p_index ++)
	{

	  //Get the coefficients for the radiative transfer:
	  
	  // 	  //1. Extinction matrix.
	  // 	  AgendaExecute{ext_mat_agenda};
	  // 	  //2. Phase Matrix.
	  // 	  AgendaExecute{pha_mat_agenda};
	  // 	  //3. Absorption Vector.
	  // 	  AgendaExecute{abs_vec_agenda};
	  // The agendas to compute the coefficients are not included yet. 
	  // For testing only
	  // dummy coefficients can be used.
	  
	  // extinction matrix
	  Matrix ext_mat(4,4,3e-12);
  
	  // absorption vector
	  Vector abs_vec(4,1e-14);
	  
	  // scattering integral vector;
	  Vector sca_vec(4,2e-12);
	  
  
	  //4. Extract sca_vec from sca_field.
	  for (Index i = 0; i < stokes_dim; i++)
	    {
	      sca_vec[i] = sca_field(p_index, 0, 0,
				      scat_za_index, 0, i);
	    }
	
	  //5. Generate Planck function.
	  Numeric T = t_field(p_index, 0, 0);
	  Numeric B;
	  Numeric f = f_grid[f_index];
	  planck(B, f, T);
	  
	  //Initialize ppath for 1D.
	  Ppath ppath_step;
	  ppath_init_structure(ppath_step, 1, 1);
  
	  // Assign value to ppath.pos:
	  ppath_step.z[0]     = z_field(p_index,0,0);
	  ppath_step.pos(0,0) = r_geoid(0,0) + ppath_step.z[0];
	  
	  // Define the direction:
	  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
	  
	  // Define the grid positions:
	  ppath_step.gp_p[0].idx   = p_index;
	  ppath_step.gp_p[0].fd[0] = 0;
	  ppath_step.gp_p[0].fd[1] = 1;
	  
	  // Calcultae intersection point of the path in the given direction
	  // with the next pressure level.
	  ppath_stepGeometric(ppath_step, 1, p_grid, lat_grid,
			      lon_grid, z_field, r_geoid, z_ground,
			      blackbody_ground);
       
  
	  // Interpolation of the intensity field on the intersection point.
	  // i_field_old corresponds to the old grid, the new grid is just
	  // one point, the intersection point.
  
	  // Calculate interpolation weights.
	  Matrix itw(1,2);
	  interpweights(itw, ppath_step.gp_p);
	  
	  // Do interpolation.
	  Vector sto_vec(4);
	  
	  // The function interp() requires Vectors as Input, Output. It can
	  // not handle numerical values. For this reason two vectors of 
	  // length 1 have to be defined.
	  Vector sto_vec_int(1);
	  Vector i_field_old_int(1);

	  // Do the interpolation for each Stokes component.
	  for (Index i = 1; i<stokes_dim; i++)
	    {
	      sto_vec_int[0] = sto_vec[i];
	      i_field_old_int[0] = i_field_old(p_index, 0, 0, 
					       scat_za_index, 0, i); 
	      interp(sto_vec_int, itw, i_field_old_int, ppath_step.gp_p);
	    }

	  // Perform RT calculation.
	  rte_scat_vecCalc(sto_vec, ext_mat, abs_vec, sca_vec, 
			   ppath_step.l_step[0], B); 
	  
	  // Assign claculated Stokes Vector to i_field. 
	  i_field(p_index, 0, 0, scat_za_index, 0, Range(joker)) = sto_vec;
  
	  // Close all loops.
	}
    }
}



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
		      const Numeric ds,
		      const Numeric B)
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




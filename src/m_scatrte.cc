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
  \file   m_scatrte.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Wed Jun 19 11:03:57 2002
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "matpackI.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "cloudbox.h"
#include "logic.h"
#include "ppath.h"
#include "agenda_class.h"
#include "physics_funcs.h"
#include "lin_alg.h"
#include "math_funcs.h"



//! Convergence test (maximum absolute difference). 
/*! 
  The function calculates the absolute differences for two successive 
  iteration fields. It picks out the maximum values for each Stokes 
  component separately. The convergence test is fullfilled under the following
  conditions:
  |I(m+1) - I(m)| <                 Intensity.
  |Q(m+1) - Q(m)| <                 The other Stokes components.
  |U(m+1) - U(m)| <  
  |V(m+1) - V(m)| <  
  These conditions have to be valid for all positions in the cloudbox and
  for all directions.
  Then convergence_flag is set to 1. 

  WS Output:
  \param convergence_flag Fag for convergence test.
  WS Input:
  \param i_field Radiation field.
  \param i_field_old Old radiation field.
  \param cloudbox_limits Limits of the cloudbox.
  \param scat_za_grid 	Zenith angle grid.
  \param scat_aa_grid Azimuth angle grid.
  \param atmosphere_dim Atmospheric dimension.

  Keyword : 
  \param epsilon   Limiting values for the convergence test
  \author Claudia Emde
  \date 2002-06-17

*/
void convergence_flagAbs(//WS Output:
                      Index& convergence_flag,
                      // WS Input:
                      const Tensor6& i_field,
                      const Tensor6& i_field_old,
                      const ArrayOfIndex& cloudbox_limits, 
                      const Vector& scat_za_grid,
                      const Vector& scat_aa_grid,
                      const Index& stokes_dim,
                      const Index& atmosphere_dim,
                      // Keyword:
                      const Vector& epsilon)
{
  assert( convergence_flag == 0 );
  
  Index N_scat_za = scat_za_grid.nelem();
  Index N_scat_aa = scat_aa_grid.nelem();

  //Loop over all components of the intensity field. 
  for (Index stokes_index = 0; stokes_index < stokes_dim; stokes_index++ )
    {
  
  if ( atmosphere_dim == 1 ) {
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
            p_index++)
        { 
          for (Index scat_za_index = 0; scat_za_index < N_scat_za;
               scat_za_index++)
            {
              Numeric diff = fabs( i_field( p_index, 0, 0, scat_za_index,
                                            0, stokes_index) -
                                   i_field_old( p_index, 0, 0, scat_za_index,
                                                0, stokes_index ));
              
              // If the absolute difference of the components is 
              // larger than the pre-defined values, return to
              // *i_fieldIterarte and do another iteration
              if( diff > epsilon[stokes_index])
                return;

            } // End loop scat_za_grid.
        }// End loop p_grid.
  } // End 1D atmosphere.
              

  //3D atmosphere:
  else if( atmosphere_dim == 3 ){
    
      for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
            p_index++)
        { 
          for (Index lat_index = cloudbox_limits[2]; lat_index <= 
                  cloudbox_limits[3]; lat_index++)
            {
              for (Index lon_index = cloudbox_limits[4]; lon_index <= 
                  cloudbox_limits[5]; lon_index++)
                {
                  for (Index scat_za_index = 0; scat_za_index < N_scat_za;
                       scat_za_index++)
                    {
                      for (Index scat_aa_index = 0; scat_aa_index < N_scat_aa;
                       scat_aa_index++)
                        {
                          Numeric diff = fabs( i_field( p_index, lat_index,
                                                        lon_index, 
                                                        scat_za_index,
                                                        scat_aa_index, 
                                                        stokes_index) -
                                               i_field_old( p_index, lat_index,
                                                            lon_index, 
                                                            scat_za_index,
                                                            scat_aa_index, 
                                                            stokes_index ));
                          
                          // If the absolute difference of the components is 
                          // larger than the pre-defined values, return to
                          // *i_fieldIterarte and do another iteration
                          if( diff > epsilon[stokes_index])
                            return;

                        }// End loop scat_aa_grid.
                    }// End loop scat_za_grid. 
                }// End loop lon_grid.
            }// End loop lat_grid. 
        }// End loop p_grid.
    
    } //end 3D atmosphere.
    
  } // end loop over Stokes dimensions. 

 // Convergence test has been successful, convergence_flag can be set to 1.
  convergence_flag = 1;
}


//! Iterative solution of the RTE.
/*! 
  A solution for the RTE with scattering is found using an iterative scheme:
  1. Calculate scattering integral.
  2. Calculate RT with fixed scattered field.
  3. Convergence test.

  WS Output:
  \param i_field       Intensity field.
  \param ppath_step    Propagation path step for RT calculation.
  \param i_field_old   Old intensity field.
  \param scat_field    Scattering integral at all grid points
                              inside the cloud box.
  \param sca_vec       Scattered field vector.
  \param stokes_vec    Stokes vector.    
  \param planck_function Planck function.
  \param l_step        Pathlength step. 
  \param convergence_flag 1 if convergence is reached after an 
                       iteration. 0 else.
  \param pha_mat_part  Scattering matrix.
  \param pha_mat_spt   Scattering matrix for a single particle type.
  \param abs_vec_spt   Absorption vector for a single particle type.
  \param ext_mat_spt   Extinction matrix for a single particle type.
                                      
  WS Input:
  \param ext_mat_agenda Agenda to compute total extinction matrix.
  \param abs_vec_agenda Agenda to compute total absorption vector.
  \param convergence_test_agenda Agenda to perform the convergence test.
  \param ppath_step_agenda Agenda to compute a propagation path step.
  \param scat_rte_agenda Agenda to compute the RTE.
  \param amp_mat       Amplitude matrix. 
  \param cloudbox_limits Limits of the cloudbox.
  \param scat_za_grid  Zenith angle grid inside the cloud box.
  \param scat_aa_grid  Azimuthal angle grid inside the cloud box.
  \param p_grid        Pressure grid.
  \param lat_grid      Latitude grid.
  \param lon_grid      Longitude grid.
  \param t_field       Temperature field for alls grid points.
  \param z_field       Geometrical altitude field.
  \param r_geoid       Matrix containing geoids.
  \param f_grid        Frequency grid.
  \param scat_f_index  Frequency index.
  \param stokes_dim    The number of Stokes components to be calculated.
  \param atmosphere_dim Atmospheric dimension.
  \param pnd_field     Particle number density field.
  \param part_types    Particle types.

  \author Claudia Emde
  \date 2002-05-29
*/
void
i_fieldIterate(
		    // WS Output:
		    Tensor6& i_field,
		    Ppath& ppath_step,
		    Tensor6& i_field_old,
		    Tensor6& scat_field,
                    Vector& sca_vec,
                    Vector& stokes_vec,
                    Numeric& planck_function,
                    Numeric& l_step,
                    Index& convergence_flag,
                    Tensor4& pha_mat_part,
                    Tensor5& pha_mat_spt,
                    Matrix& abs_vec_spt,
                    Tensor3& ext_mat_spt,
                    // WS Input:
                    const Agenda& ext_mat_agenda,
                    const Agenda& abs_vec_agenda,
                    const Agenda& convergence_test_agenda,
                    const Agenda& ppath_step_agenda,
                    const Agenda& scat_rte_agenda,
		    const Tensor6& amp_mat,
		    const ArrayOfIndex& cloudbox_limits,
		    const Vector& scat_za_grid,
		    const Vector& scat_aa_grid,
		    const Vector& p_grid,
		    const Vector& lat_grid,
		    const Vector& lon_grid,
		    const Tensor3& t_field,
		    const Tensor3& z_field,
		    const Matrix& r_geoid,
		    const Vector& f_grid,
		    const Index& scat_f_index,
                    const Index& stokes_dim,
                    const Index& atmosphere_dim,
                    const Tensor4& pnd_field,
                    const Vector& part_types
                    )
{
  // Check the input
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
			"The dimension of stokes vector must be"
			"1,2,3, or 4");
 
  if (atmosphere_dim == 3){
 
   //Does i_field have the right dimension? 
    assert( is_size( i_field, p_grid.nelem(), lat_grid.nelem(), 
		      lon_grid.nelem(), scat_za_grid.nelem(), 
                     scat_aa_grid.nelem(), stokes_dim));
    
    //Dimension of cloudbox_limits
    assert(is_size(cloudbox_limits, 6));
    }

  else if (atmosphere_dim == 1 ){
    assert ( is_size( i_field, p_grid.nelem(), 1, 
		      1, scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));
  }
  
  else if (atmosphere_dim == 2){
    throw runtime_error(
                        "Scattering calculations are not possible for a 2D"
                        "atmosphere. If you want to do scattering calculations"
                        "*atmosphere_dim* has to be either 1 or 3"
                        );
      }

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  // Grids have to be adapted to atmosphere_dim.
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Is the frequency index valid?
  assert( scat_f_index <= f_grid.nelem() );

  
  
  //The following steps are repeated until convergence is reached.
  convergence_flag = 0;
  while(convergence_flag == 0) {
  
    // 1. Copy i_field to i_field_old.
    i_field_old = i_field;
    
    
    // 2.Calculate scattered field vector for all points in the cloudbox.
    
    // Resize variables:
    ext_mat_spt.resize(part_types.nelem(), stokes_dim, stokes_dim);
    abs_vec_spt.resize(part_types.nelem(), stokes_dim);
    pha_mat_spt.resize(part_types.nelem(), scat_za_grid.nelem(), 1,
                       stokes_dim, stokes_dim);
    
    // To calculate the scattering integral, pha_mat_spt is required 
    // for all directions. Furthermore is is required to calculate 
    // the absortion vector.
    
    Index N_scat_za = scat_za_grid.nelem();
    for(Index scat_za_index = 0; scat_za_index < N_scat_za;
        scat_za_index ++)
      {
        ext_mat_sptCalc(ext_mat_spt, amp_mat, scat_za_index, 0, 
                        scat_f_index, f_grid);
        pha_mat_sptCalc(pha_mat_spt, amp_mat, scat_za_index, 0);
      }
    
    // Total scattering matrix.
    //pha_mat_partCalc(pha_mat_part, pha_mat_spt, pnd_field, 
    //cloudbox_limits, atmosphere_dim);
    
    
    scat_field = i_field;
    // ---- here will be the function to calculate the scattering integral
  
    
    //Update i_field.
    if( atmosphere_dim == 1 )
      {
        i_fieldUpdate1D(//Output:
                        i_field, ppath_step, stokes_vec, 
                        sca_vec, planck_function, l_step,  
                        abs_vec_spt, ext_mat_spt, 
                        //Input:
                        ext_mat_agenda, abs_vec_agenda, ppath_step_agenda,
                        scat_rte_agenda,  amp_mat, scat_field,
                        cloudbox_limits, scat_za_grid, p_grid, 
                        t_field, z_field, r_geoid, f_grid, scat_f_index,  
                        pnd_field, stokes_dim, atmosphere_dim, part_types,
                        pha_mat_spt);
      }
    
    //Convergence test.
    convergence_test_agenda.execute();
  }//end of while loop, convergence is reached.
}


//! 1D RT calculation inside the cloud box.
/*! 
  This function loops over all grid points and all directions and performs 
  the RT calculation with a fixed scattering integral for one frequency 
  of the frequency grid specified by *scat_f_index*. 
  
  The loop over directions is the outermost loop. Here the optical properties
  for single particle types are calculated as they are not depending on the 
  position of the particles. 
  The inner loop is the loop over the positions. Inside this loop the total 
  optical properties including the partile types as well as the gaseous
  species are calculated. Then the radiative transfer equation can be computed.

  WS Output:
  \param i_field       Updated intensity field. 
  \param ppath_step    Propagation path step for RT calculation.
  \param stokes_vec    Stokes vector.
  \param sca_vec       Scattered field vector.
  \param planck_function Planck function.
  \param l_step        Pathlength step.
  \param abs_vec_spt   Absorption vector for a single particle type.
  \param ext_mat_spt   Extinction matrix for a single particle type.
  WS Input:
  \param ext_mat_agenda Agenda to calculate extinction matrix.
  \param abs_vec_agenda Agenda to calculate absorption vector.
  \param ppath_step_agenda Agenda to compute a propagation path step.
  \param scat_rte_agenda Agenda to compute the RTE.
  \param amp_mat       Amplitude matrix. 
  \param scat_field     Scattering integral at all grid points
                              inside the cloud box.
  \param cloudbox_limits Limits of the cloudbox.
  \param scat_za_grid  Zenith angle grid inside the cloud box.
  \param p_grid        Pressure grid (is required only for checking the 
                       input).
  \param t_field       Temperature field for all grid points.
  \param z_field       Geometrical altitude field.
  \param r_geoid       Matrix containing geoids.
  \param f_grid        Frequency grid.
  \param scat_f_index  Frequency index.
  \param pnd_field     Particle number density field.
  \param stokes_dim    The number of Stokes components to be calculated.
  \param atmosphere_dim Dimension of the atmosphere.
  \param part_types    Particle types.
  \param pha_mat_spt   Scattering matrix for a single particle type.
  

Some workspace variables have to be defined to execute the agendas. These are:
         ext_mat       Extinction matrix (4x4 matrix).
         abs_vec       Absorprion vector (4 elements).
         ext_mat_part  Extinction matrix only for particles 
                       (no gaseous extinction).
         abs_vec_part  Absorprion vector only for particles.             
         blackbody_ground Needed for the ppath_step calculation.  

  \author Claudia Emde
  \date 2002-05-30
*/
void
i_fieldUpdate1D(// WS Output:
		Tensor6& i_field,
		Ppath& ppath_step,
                Vector& stokes_vec,
                Vector& sca_vec,
                Numeric& planck_function,
                Numeric& l_step,
                Matrix& abs_vec_spt,
                Tensor3& ext_mat_spt,
                // WS Input:
                const Agenda& ext_mat_agenda,
                const Agenda& abs_vec_agenda,
		const Agenda& ppath_step_agenda,
                const Agenda& scat_rte_agenda,
                const Tensor6& amp_mat,
		const Tensor6& scat_field,
		const ArrayOfIndex& cloudbox_limits,
		const Vector& scat_za_grid,
                const Vector& p_grid,
                const Tensor3& t_field,
		const Tensor3& z_field,
		const Matrix& r_geoid,
		const Vector& f_grid,
		const Index& scat_f_index,
		const Tensor4& pnd_field,
                const Index& stokes_dim,
                const Index& atmosphere_dim,
                const Vector& part_types,
                const Tensor5& pha_mat_spt
                )
{

  //Check the input
  
  assert( atmosphere_dim == 1);
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert ( is_size( i_field, p_grid.nelem(), 1, 
                    1, scat_za_grid.nelem(), 
                    1, stokes_dim));
  
  assert ( is_size( scat_field , p_grid.nelem(), 1, 
                    1, scat_za_grid.nelem(), 
                    1, stokes_dim));  


  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();


 
  //Loop over all directions, defined by scat_za_grid 
  for(Index scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {

      //Calculate optical properties for single particle types:

      //Calculate ext_mat_spt for the direction 
      //corresponding to the outer loop:
      ext_mat_sptCalc(ext_mat_spt, amp_mat, scat_za_index, 0, scat_f_index, 
                      f_grid);      

      // For a 1D atmosphere the azimuthal angle grid is not defined. 
      // Only 1 value, which is arbitrary set to 0, is passed into the function
      // abs_vec_sptCalc. 
      // pha_mat_spt is already calculated in i_fieldIterate.
      Vector scat_aa_grid(1);
      scat_aa_grid[0] = 0;
      abs_vec_sptCalc(abs_vec_spt, ext_mat_spt, pha_mat_spt, scat_za_grid, 
                      scat_aa_grid);
    

      //Loop over all positions inside the cloud box defined by scat_p_grid
      for(Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	  p_index ++)
	{
	  //Print the loop indices (just for testing the function)

	  cout << "\n loop indices: \n";
	  cout << "\n scat_za_index ---------"<< scat_za_index;
	  cout << "\n p_index       ---------"<< p_index;
          cout << "\n stokes_dim    ---------"<< stokes_dim;

	  //Get the coefficients for the radiative transfer:
          
          //1.Extinction Matrix.
          ext_mat_agenda.execute();

          //3. Absorption Vector.
	  abs_vec_agenda.execute();


          // Get sca_vec and stokes_vec from the fields.
          sca_vec.resize(stokes_dim);
          stokes_vec.resize(stokes_dim);
          
	  for (Index i = 0; i < stokes_dim; i++)
	    {
              //Extract sca_vec from sca_field.
	      sca_vec[i] = scat_field(p_index, 0, 0,
				      scat_za_index, 0, i);
              //Extract stokes_vec from i_field.
              stokes_vec[i] = i_field(p_index, 0, 0,
                                         scat_za_index, 0, i);
	     }
	

	  //Generate Planck function.
	  Numeric T1 = t_field(p_index, 0, 0);
          Numeric T2 = t_field(p_index + 1, 0, 0); 
	  Numeric f = f_grid[scat_f_index];
	  planck_function = 0.5*(planck(f, T1)+planck(f,T2));

	  //Initialize ppath for 1D.
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
	  
	  // Call ppath_step_agenda: 
          ppath_step_agenda.execute();
	  
	  // Check if the agenda has returned ppath.step with reasonable 
	  // values. 
          //	  cout << "\n ";
	  //PpathPrint( ppath_step, "ppath");
          
          l_step = ppath_step.l_step[0];
          
          // Call scat_rte_agenda:
          scat_rte_agenda.execute();

	  // Assign calculated Stokes Vector to i_field. 
	  i_field(p_index, 0, 0, scat_za_index, 0, Range(joker)) = stokes_vec;
  
	  // Close all loops.
	}
    }
}



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
  \param planck_function  Input: Planck function.
  \param stokes_dim Input: Stokes dimension. 

  \author Claudia Emde
  \date 2002-06-08
*/
void
stokes_vecGeneral(//WS Output and Input:
               Vector& stokes_vec,
               //WS Input:
               const Matrix& ext_mat,
               const Vector& abs_vec,
               const Vector& sca_vec,
               const Numeric& l_step,
               const Numeric& planck_function,
               const Index& stokes_dim)
{ 
  // Stokes dimension
  assert(stokes_dim <= 4 && stokes_dim != 0);
 
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
  B_abs_vec *= planck_function; 
  
  for (Index i=0; i<stokes_dim; i++) 
    b[i] = abs_vec[i] + sca_vec[i];  // b = abs_vec * B + sca_vec

  // solve K^(-1)*b = x
  ludcmp(LU, indx, ext_mat);
  lubacksub(x, LU, b, indx);

  Matrix ext_mat_ds(stokes_dim, stokes_dim);
  ext_mat_ds = ext_mat;
  ext_mat_ds *= -l_step; // ext_mat = -ext_mat*ds

  Index q = 10;  // index for the precision of the matrix exponential function
  Matrix exp_ext_mat(stokes_dim, stokes_dim);
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
  \param planck_function  Input: Planck function.
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
	      const Numeric& planck_function,
	      const Index& stokes_dim)
{ 
  //Check if we really consider the scalar case.
  if(stokes_dim != 1)
    throw runtime_error("You can use the method *stokes_vecScalar* only if"
                        "*stokes_dim* equals 1"); 

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
    (abs_coeff*planck_function+sca_int1D) /
    ext_coeff* (1-exp(-ext_coeff*l_step));
 
  //Put the first component back into *sto_vec*:
  stokes_vec[0] = stokes_vec1D;
}
		


//! This method computes the scattering integral

/*! 
By scattering integral, we mean the field generated by integrating
the product of intensity field and phase matrix over all incident 
angles.  
  
\param sca_vec Output :scattering integraal
\param i_field Input : the intensity field 
\param pha_mat_part Input : the phase matrix
\param scat_za_grid zenith angle grid
\param scat_aa_grid azimuth angle grid
\param p_grid pressure grid
\param lat_grid latitude grid
\param lon_grid longitude grid
\param stokes_dim dimension of stokes vector
\param atmosphere_dim atmospheric dimension
*/
void
scat_integralCalc(//WS Output:
              Tensor6& scat_field,
              //WS Input: 
	      const Tensor6& i_field,
	      const Tensor4& pha_mat_part,
	      const Vector& scat_za_grid,
	      const Vector& scat_aa_grid,
	      const Vector& p_grid,
	      const Vector& lat_grid,
	      const Vector& lon_grid,
	      const Index& stokes_dim,
	      const Index& atmosphere_dim,
	      const ArrayOfIndex& cloudbox_limits
	      )
  
{
  Index Nza = scat_za_grid.nelem();
  Index Naa = scat_aa_grid.nelem();

  Matrix product_field(stokes_dim, stokes_dim);
  if (atmosphere_dim ==3){
    
    assert ( is_size( i_field, p_grid.nelem(), lat_grid.nelem(), 
		      lon_grid.nelem(), scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));
    assert ( is_size( scat_field, p_grid.nelem(), lat_grid.nelem(), 
		      lon_grid.nelem(), scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));  
  }
  
  else if (atmosphere_dim == 1 ){
    assert ( is_size( i_field, p_grid.nelem(), 1, 
		      1, scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));
    assert ( is_size( scat_field, p_grid.nelem(), 1, 
		      1, scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));
  }
    
  if( atmosphere_dim == 1 ){
    
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++)
      {
	for (Index za_index = 0; za_index < Nza; ++ za_index)
	  {
	    for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
	      {
		for (Index i = 0; i < stokes_dim; i++)
		  {
		    ConstMatrixView pha = pha_mat_part(za_index,
						       aa_index,
						       Range(joker),
						       Range(joker));
		    ConstVectorView ifield_in = i_field(p_index,
							1,1,
							za_index,
							aa_index,
							Range(joker));

		      mult(product_field, pha, ifield_in);
		    

		    scat_field
		      = AngIntegrate_trapezoid(product_field,
					       scat_za_grid,
					       scat_aa_grid);
		    
		  }
	      }
	  }
      }
  }
}

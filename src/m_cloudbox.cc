/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                      Claudia Emde     <claudia@sat.physik.uni-bremen.de>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_cloudbox.cc
  \author Patrick Eriksson and Claudia Emde
  \date   2002-05-08 

  \brief  Workspace functions related to the cloud box.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.*/



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

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! CloudboxOff
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void CloudboxOff(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits )
{
  cloudbox_on = 0;
  cloudbox_limits.resize(0);
}



//! CloudboxSetManually
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-19
*/
void CloudboxSetManually(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
        // WS Input:
        const Index&    atmosphere_dim,
        const Vector&   p_grid,
        const Vector&   lat_grid,
        const Vector&   lon_grid,
        const Index&    blackbody_ground,
        // Control Parameters:
        const Numeric& p1,
        const Numeric& p2,
        const Numeric& lat1,
        const Numeric& lat2,
        const Numeric& lon1,
        const Numeric& lon2 )
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Check keyword arguments
  if( p1 <= p2 )
    throw runtime_error( 
            "The pressure in *p1* must be bigger than the pressure in *p2*." );
  if( p1 <= p_grid[p_grid.nelem()-1] )
    throw runtime_error( "The pressure in *p1* must be larger than the "
                                                   "last value in *p_grid*." );
  if( p2 >= p_grid[0] )
    throw runtime_error( "The pressure in *p2* must be smaller than the "
                                                  "first value in *p_grid*." );
  if( atmosphere_dim >= 2 )
    {
      if( lat2 <= lat1 )
	throw runtime_error( 
	 "The latitude in *lat2* must be bigger than the latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
	throw runtime_error( "The latitude in *lat1* must be >= than the "
                                               "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
	throw runtime_error( "The latitude in *lat2* must be <= than the "
                                         "next to last value in *lat_grid*." );
    }
  if( atmosphere_dim == 3 )
    {
      if( lon2 <= lon1 )
	throw runtime_error( 
       "The longitude in *lon2* must be bigger than the longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
	throw runtime_error( "The longitude in *lon1* must be >= than the "
                                               "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
	throw runtime_error( "The longitude in *lon2* must be <= than the "
                                         "next to last value in *lon_grid*." );
    }

  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize( atmosphere_dim*2 );

  // Pressure limits
  if( p1 > p_grid[1] )
    {
      cloudbox_limits[0] = 0;
    }
  else
    {
      for( cloudbox_limits[0]=1; p_grid[cloudbox_limits[0]+1]>p1; 
                                                     cloudbox_limits[0]++ ) {}
    }
  if( !blackbody_ground && cloudbox_limits[0]!=0 )
    {
      ostringstream os;
      os << "The lower vertical limit of the cloud box must be the lowest "
         << "pressure\nsurface when the ground is not a blackbody.";
      throw runtime_error( os.str() );
    }
  if( p2 < p_grid[p_grid.nelem()-2] )
    {
      cloudbox_limits[1] = p_grid.nelem() - 1;
    }
  else
    {
      for( cloudbox_limits[1]=p_grid.nelem()-2; 
                    p_grid[cloudbox_limits[1]-1]<p2; cloudbox_limits[1]-- ) {}
    }

  // Latitude limits
  if( atmosphere_dim >= 2 )
    {
      for( cloudbox_limits[2]=1; lat_grid[cloudbox_limits[2]+1]<lat1; 
                                                     cloudbox_limits[2]++ ) {}
      for( cloudbox_limits[3]=lat_grid.nelem()-2; 
                lat_grid[cloudbox_limits[3]-1]>lat2; cloudbox_limits[3]-- ) {}
    }

  // Longitude limits
  if( atmosphere_dim == 3 )
    {
      for( cloudbox_limits[4]=1; lon_grid[cloudbox_limits[4]+1]<lon1; 
                                                     cloudbox_limits[4]++ ) {}
      for( cloudbox_limits[5]=lon_grid.nelem()-2; 
                lon_grid[cloudbox_limits[5]-1]>lon2; cloudbox_limits[5]-- ) {}
    }
}



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

  \author Claudia Emde
  \date 2002-06-17

*/
void convergence_flagAbs(//WS Output:
                      Index& convergence_flag,
                      //WS Input:
                      const Tensor6& i_field,
                      const Tensor6& i_field_old,
                      const ArrayOfIndex& cloudbox_limits, 
                      const Vector& scat_za_grid,
                      const Vector& scat_aa_grid,
                      const Index& stokes_dim)
{
  assert( convergence_flag == 0 );
  
  Index N_scat_za = scat_za_grid.nelem();
  Index N_scat_aa = scat_aa_grid.nelem();

  // define the limiting values
  Vector epsilon(4);
  epsilon[0] = 1e-20;
  epsilon[1] = 1e-21;
  epsilon[2] = 1e-21;
  epsilon[3] = 1e-21;

  //Loop over all components of the intensity field. 
  for (Index stokes_index = 0; stokes_index < stokes_dim; stokes_index++ )
    {
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
                                                        lon_index, scat_za_index,
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
                        }
                    }
                }
            }
        }
    }

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
  \param scat_field     Scattering integral at all grid points
                              inside the cloud box.
  \param ext_mat       Extinction matrix.
  \param abs_vec       Absorption vector.
  \param sca_vec       Scattered field vector.
  \param stokes_vec    Stokes vector.    
  \param planck_function Planck function.
  \param l_step        Pathlength step. 
  \param convergence_flag 1 if convergence is reached after an 
                       iteration. 0 else.
                            
  WS Input:
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
  \param z_ground      Ground altitude.
  \param r_geoid       Matrix containing geoids.
  \param f_grid        Frequency grid.
  \param scat_f_index  Frequency index.
  \param part_types    Vector containing particle types.
  \param blackbody_ground Flag to treat ground as blackbody.
  \param stokes_dim    The number of Stokes components to be calculated.
  \param atmosphere_dim Atmospheric dimension.

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
                    Matrix& ext_mat,
                    Vector& abs_vec,
                    Vector& sca_vec,
                    Vector& stokes_vec,
                    Numeric& planck_function,
                    Numeric& l_step,
                    Index& convergence_flag,
		    // WS Input:
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
		    const Matrix& z_ground,
		    const Matrix& r_geoid,
		    const Vector& f_grid,
		    const Index& scat_f_index,
		    const Vector& part_types,
		    const Index& blackbody_ground,
		    const Index& stokes_dim,
                    const Index& atmosphere_dim
                    )
{
  // Check the input
  
  assert( is_size( amp_mat, part_types.nelem(), scat_za_grid.nelem(), 
		   scat_aa_grid.nelem(), scat_za_grid.nelem(),
		   scat_aa_grid.nelem(), 8));
  
  if (atmosphere_dim == 3){
    assert ( is_size( i_field, p_grid.nelem(), lat_grid.nelem(), 
		      lon_grid.nelem(), scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));
  }
  else if (atmosphere_dim == 1 ){
    assert ( is_size( i_field, p_grid.nelem(), 1, 
		      1, scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(), stokes_dim));
  }

  //The following steps are repeated until convergence is reached.
  convergence_flag = 0;
  while(convergence_flag == 0) {
  

  // Copy i_field to i_field_old.
   i_field_old = i_field;
  
  //Calculate scattered field vector for all points in the cloudbox.

  scat_field = i_field;
  // ---- here will be the function to calculate the scattering integral
  
  
  //Update i_field.
  if( atmosphere_dim == 1 )
    {
      i_fieldUpdate1D(//Output:
                      i_field, ppath_step, stokes_vec, ext_mat, abs_vec, sca_vec,
                      planck_function, l_step,
                      //Input:
                      ppath_step_agenda, scat_rte_agenda,
                      i_field_old, amp_mat, scat_field,
		      cloudbox_limits, scat_za_grid, scat_aa_grid, p_grid, 
		      lat_grid, lon_grid, t_field, z_field, z_ground,
		      r_geoid, f_grid, scat_f_index, blackbody_ground, 
		      stokes_dim);
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

  Note: The function uses the same input and output variables as the
  equivalent function for the 
  3D calculations. Here the dimensions which are not needed are empty.
  
  The agendas for computing ext_mat, pha_mat and abs_vec are needed as further
  input!
  
  WS Output:
  \param i_field       Updated intensity field. 
  \param ppath_step    Propagation path step for RT calculation.
  \param stokes_vec    Stokes vector.
  \param ext_mat       Extinction matrix.
  \param abs_vec       Absorption vector.
  \param sca_vec       Scattered field vector.
  \param planck_function Planck function.
  \param l_step        Pathlength step.
  WS Input:
  \param ppath_step_agenda Agenda to compute a propagation path step.
  \param scat_rte_agenda Agenda to compute the RTE.
  \param i_field_old   Old intensity field.
  \param amp_mat       Amplitude matrix. 
  \param scat_field     Scattering integral at all grid points
                              inside the cloud box.
  \param cloudbox_limits Limits of the cloudbox.
  \param scat_za_grid  Zenith angle grid inside the cloud box.
  \param scat_aa_grid  Azimuthal angle grid inside the cloud box.
  \param p_grid        Pressure grid.
  \param lat_grid      Latitude grid.
  \param lon_grid      Longitude grid.
  \param t_field       Temperature field for all grid points.
  \param z_field       Geometrical altitude field.
  \param z_ground      Ground altitude.
  \param r_geoid       Matrix containing geoids.
  \param f_grid        Frequency grid.
  \param scat_f_index  Frequency index.
  \param blackbody_ground Flag to treat ground as blackbody.
  \param stokes_dim    The number of Stokes components to be calculated.

  \author Claudia Emde
  \date 2002-05-30
*/
void
i_fieldUpdate1D(// WS Output:
		Tensor6& i_field,
		Ppath& ppath_step,
                Vector& stokes_vec,
                Matrix& ext_mat,
                Vector& abs_vec,
                Vector& sca_vec,
                Numeric& planck_function,
                Numeric& l_step,
		// WS Input:
		const Agenda& ppath_step_agenda,
                const Agenda& scat_rte_agenda,
		const Tensor6& i_field_old,
		const Tensor6& amp_mat,
		const Tensor6& scat_field,
		const ArrayOfIndex& cloudbox_limits,
		const Vector& scat_za_grid,
		const Vector& scat_aa_grid,
		const Vector& p_grid,
		const Vector& lat_grid,
		const Vector& lon_grid,
		const Tensor3& t_field,
		const Tensor3& z_field,
		const Matrix& z_ground,
		const Matrix& r_geoid,
		const Vector& f_grid,
		const Index& scat_f_index,
		const Index& blackbody_ground,
		const Index& stokes_dim
             	)
{

  //Check the input
 assert ( is_size( i_field, p_grid.nelem(), 1, 
		   1, scat_za_grid.nelem(), 
		   scat_aa_grid.nelem(), stokes_dim));
  
 assert ( is_size( i_field_old, p_grid.nelem(), 1, 
		   1, scat_za_grid.nelem(), 
		   scat_aa_grid.nelem(), stokes_dim));  
  
 assert ( is_size( scat_field , p_grid.nelem(), 1, 
		   1, scat_za_grid.nelem(), 
		   scat_aa_grid.nelem(), stokes_dim));  


  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
 
  //Loop over all directions, defined by scat_za_grid 
  for(Index scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
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
	  
	  // 	  //1. Extinction matrix.
	  // 	  AgendaExecute{ext_mat_agenda};
	  // 	  //2. Phase Matrix.
	  // 	  AgendaExecute{pha_mat_agenda};
	  // 	  //3. Absorption Vector.
	  // 	  AgendaExecute{abs_vec_agenda};
	  // The agendas to compute the coefficients are not included yet. 
	  // For testing only
	  // dummy coefficients can be used.
	  
	  //extinction matrix
	  Matrix A(4,4,3e-12);
          ext_mat.resize(stokes_dim, stokes_dim);
          ext_mat = A;
          
          //absorption vector
	  Vector b(4,1e-14);
          abs_vec.resize(stokes_dim); 
	  abs_vec = b;
         
          //scattering integral vector;
	  Vector s(4,1);
          sca_vec.resize(stokes_dim);
          stokes_vec.resize(stokes_dim);

          stokes_vec = s;
         
	  //Extract sca_vec from sca_field.
	  for (Index i = 0; i < stokes_dim; i++)
	    {
	      sca_vec[i] = scat_field(p_index, 0, 0,
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
	  cout << "\n ";
	  PpathPrint( ppath_step, "ppath");
          
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
  Numeric stokes_vec1D = stokes_vec[0];
  
  //Do a radiative transfer step calculation:
  stokes_vec1D = stokes_vec1D * exp(-ext_coeff*l_step) + 
    (abs_coeff*planck_function+sca_int1D) /
    ext_coeff* (1-exp(-ext_coeff*l_step));
 
  //Put the first component back into *sto_vec*:
  stokes_vec[0] = stokes_vec1D;
}
		



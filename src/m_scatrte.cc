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
  \author Sreerekha T.R. <rekha@uni-bremen.de>
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
#include <cmath>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "cloudbox.h"
#include "logic.h"
#include "ppath.h"
#include "agenda_class.h"
#include "physics_funcs.h"
#include "lin_alg.h"
#include "math_funcs.h"
#include "messages.h"
#include "xml_io.h"

extern const Numeric PI;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


//! Convergence test (maximum absolute difference). 
/*! 
  The function calculates the absolute differences for two successive 
  iteration fields. It picks out the maximum values for each Stokes 
  component separately. The convergence test is fullfilled under the following
  conditions:
  |I(m+1) - I(m)| < epsilon_1                Intensity.
  |Q(m+1) - Q(m)| < epsilon_2                The other Stokes components.
  |U(m+1) - U(m)| < epsilon_3 
  |V(m+1) - V(m)| < epsilon_4
  The limits for convergence have to be set in the controlfile by setting the 
  vector *epsilon* to appropriate values.

  The conditions have to be valid for all positions in the cloudbox and
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
 
  //Check the input:
  assert( convergence_flag == 0 );
  assert( atmosphere_dim == 1 || atmosphere_dim == 3 );
  assert( cloudbox_limits.nelem() == 2*atmosphere_dim);
  assert( 0 < stokes_dim < 5 );

  // Check keyword "epsilon":
  if ( epsilon.nelem() != stokes_dim )
    throw runtime_error(
                        "You have to specify limiting values for the "
                        "convergence test for each Stokes component "
                        "separately. That means that *epsilon* must "
                        "have *stokes_dim* elements!"
                        );
  
  static Index counter = 0;
  counter = counter+1;
  cout << "Number of iterations:   "<< counter << endl;
  
  ostringstream os;

  os << counter;
  xml_write_to_file("i_field_" + os.str() + ".xml", i_field);

  cout<< "i_field    " << i_field << endl;

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
              Numeric diff = 
		fabs( i_field(( p_index-cloudbox_limits[0]), 
			      0, 0, scat_za_index, 0, stokes_index) -
		      i_field_old(( p_index-cloudbox_limits[0]),
				  0, 0, scat_za_index, 0, stokes_index ));
             
	      
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
                          Numeric diff =
			    fabs( i_field((p_index-cloudbox_limits[0]),
					  (lat_index-cloudbox_limits[2]),
					  (lon_index-cloudbox_limits[4]), 
					  scat_za_index,
					  scat_aa_index, 
					  stokes_index) -
				  i_field_old((p_index-cloudbox_limits[0]),
					      (lat_index-cloudbox_limits[2]),
					      (lon_index-cloudbox_limits[4]), 
					      scat_za_index,
					      scat_aa_index, 
					      stokes_index ));
                          
                          // If the absolute difference of the components is 
                          // larger than the pre-defined values, return to
                          // *i_fieldIterarte* and do another iteration
			  
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
       To calculate the scattering integral the total scattering matrix is 
       required. This is obtained using the functions *pha_mat_sptCalc* and 
       *pha_matCalc*. The method *scat_fieldCalc* performs the integration.
  2. Calculate RT with fixed scattered field.
       The radiative transfer equation with fixed scattering integral can be 
       solved analytically if the coefficients are assumed to be constant.
       According to *atmosphere_dim* either *i_fieldUpdate1D* or 
       *i_fieldUpdate3D* are called to perform the calculation. Inside these
       methods the agenda *scat_rte_agenda* is executed.  
  3. Convergence test.
       Here the *convergence_test_agenda* is executed.

  Note: The atmospheric dimensionality *atmosphere_dim* can be either 1 or 3. 
        To these dimensions the method adapts automatically. 
        If *atmosphere_dim* equals 2, it returns an error message, as 2D
        scattering calculations can not be performed.

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
  \param pha_mat  Scattering matrix.
  \param pha_mat_spt   Scattering matrix for a single particle type.
  \param abs_vec_spt   Absorption vector for a single particle type.
  \param ext_mat_spt   Extinction matrix for a single particle type.
  \param scat_p_index  Pressure index.
  \param scat_lat_index Latitude index.
  \param scat_lon_index Longitude index.
                                      
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
                    Tensor4& pha_mat,
                    Tensor5& pha_mat_spt,
                    Matrix& abs_vec_spt,
                    Tensor3& ext_mat_spt,
                    Matrix& ext_mat,
                    Vector& abs_vec,
                    Index& scat_p_index,
                    Index& scat_lat_index,
                    Index& scat_lon_index,
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
    
    // Does i_field have the right dimension? 
    assert ( is_size( i_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
		      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
		      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
		      scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(),
		      stokes_dim));
    
    //Dimension of cloudbox_limits
    assert(is_size(cloudbox_limits, 6));
  }
 
  else if (atmosphere_dim == 1 ){
    assert ( is_size( i_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
		      1, 
		      1,
		      scat_za_grid.nelem(), 
		      1,
		      stokes_dim));
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
    pha_mat_spt.resize(part_types.nelem(), scat_za_grid.nelem(), 
                       scat_aa_grid.nelem(), stokes_dim, stokes_dim);
    
    // To calculate the scattering integral, pha_mat_spt is required 
    // for all directions. Furthermore is is required to calculate 
    // the absortion vector.
   
    Index N_scat_za = scat_za_grid.nelem();
    Index N_scat_aa = scat_aa_grid.nelem();
    
    for(Index scat_za_index = 0; scat_za_index < N_scat_za;
        scat_za_index ++)
      {
        //If we have i_field_dim == 1, then we can avoid this loop, 
        //it would always contribute with a factor 2*pi (CE). 
        //This should be discussed ...
        
	//The loop over azimuth angle (STR).
	for(Index scat_aa_index = 0; scat_aa_index < N_scat_aa;
	    scat_aa_index ++)
	  {
	    ext_mat_sptCalc(ext_mat_spt, 
			    amp_mat,
			    scat_za_index,
			    scat_aa_index,
			    scat_f_index,
			    f_grid);
	    
	    pha_mat_sptCalc(pha_mat_spt,
			    amp_mat,
			    scat_za_index,
			    scat_aa_index);
	  }
      }
    
    // Calculate the scattered field.
    scat_fieldCalc(scat_field, pha_mat, i_field, pha_mat_spt, pnd_field, 
                   scat_za_grid, scat_aa_grid, p_grid, lat_grid, lon_grid,
                   stokes_dim, atmosphere_dim, cloudbox_limits);
    
    // Update i_field.
    if( atmosphere_dim == 1 )
      {
        i_fieldUpdate1D(//Output:
                        i_field, ppath_step, stokes_vec, 
                        sca_vec, planck_function, l_step,  
                        abs_vec_spt, ext_mat_spt,ext_mat, abs_vec,
                        scat_p_index,
                        //Input:
                        ext_mat_agenda, abs_vec_agenda, ppath_step_agenda,
                        scat_rte_agenda,  amp_mat, scat_field,
                        cloudbox_limits, scat_za_grid, scat_aa_grid, p_grid, 
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
  \param ext_mat       Extinction matrix (4x4 matrix).
  \param abs_vec       Absorprion vector (4 elements).
  \param scat_p_index  Pressure index.
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
  \param scat_aa_grid  Azimuth angle grid inside the cloud box.//STR
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
                Matrix& ext_mat,
                Vector& abs_vec,
                Index& scat_p_index,
                // WS Input:
                const Agenda& ext_mat_agenda,
                const Agenda& abs_vec_agenda,
		const Agenda& ppath_step_agenda,
                const Agenda& scat_rte_agenda,
                const Tensor6& amp_mat,
		const Tensor6& scat_field,
		const ArrayOfIndex& cloudbox_limits,
		const Vector& scat_za_grid,
		const Vector& scat_aa_grid,//STR
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
  
  assert ( is_size( i_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
		      1, 
		      1,
		      scat_za_grid.nelem(), 
		      1,
		      stokes_dim));
  assert ( is_size( scat_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
		      1, 
		      1,
		      scat_za_grid.nelem(), 
		      1,
		      stokes_dim));  
  
  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  const Index N_scat_aa = scat_aa_grid.nelem();

  //Loop over all directions, defined by scat_za_grid 
  for(Index scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      //Loop over all directions, defined by scat_aa_grid 
      for(Index scat_aa_index = 0; scat_aa_index < N_scat_aa; scat_aa_index ++)
	{
	
      //Calculate optical properties for single particle types:
    
      //Calculate ext_mat_spt for the direction 
      //corresponding to the outer loop:
      ext_mat_sptCalc(ext_mat_spt,
		      amp_mat,
		      scat_za_index,
		      scat_aa_index,
		      scat_f_index, 
                      f_grid);   
   
      // cout << "ext_mat_spt is calculated"<<"\n";
      // For a 1D atmosphere the azimuthal angle grid is not defined. 
      // Only 1 value, which is arbitrary set to 0, is passed into the function
      // abs_vec_sptCalc. 
      // pha_mat_spt is already calculated in i_fieldIterate.
      // next two line commented by STR
      // Vector scat_aa_grid(1);
      // scat_aa_grid[0] = 0;

      // (CE): If the radiation field and the atmodpheric profiles are 1D
      //       the calculations will be faster if we use a 1D intergation
      //       routine as the integration over the azimutal angle gives 
      //       always 2*pi.

      abs_vec_sptCalc(abs_vec_spt, 
		      ext_mat_spt,
		      pha_mat_spt,
		      scat_za_grid, 
                      scat_aa_grid);
     
      //cout << "abs_vec_spt is calculated"<<"\n";
	} // close loop over scat_aa_grid
    }// close loop over scat_za_grid
   
  

  for(Index scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    { 
      
      // Calculate ext_mat, abs_vec for all points inside the cloudbox.
      // sca_vec can be obtained from the workspace variable scat_field.
      // As we need the average for each layer, it makes sense to calculte
      // the coefficients once and store them in an array instead of 
      // calculating at each point the coefficient of the point above and 
      // the point below. 
     
      ArrayOfMatrix ext_mat_array;
      ext_mat_array.resize(cloudbox_limits[1]-cloudbox_limits[0]+1);

      ArrayOfVector abs_vec_array;
      abs_vec_array.resize(cloudbox_limits[1]-cloudbox_limits[0]+1);

      Vector planck_vector;
      planck_vector.resize(cloudbox_limits[1]-cloudbox_limits[0]+1);

     
      // Loop over all positions inside the cloudbox defined by the 
      // cloudbox_limits.
      for(Index p_index = cloudbox_limits[0]; p_index
	<= cloudbox_limits[1]; p_index ++)
	{
	  
	  //Print the loop indices (just for testing the function)
	  
	  cout << "\n loop indices: \n";
	  cout << "\n scat_za_index ---------"<< scat_za_index;
	  cout << "\n p_index       ---------"<< p_index;
          cout << "\n stokes_dim    ---------"<< stokes_dim;
	  cout << "\n cloudbox_limits    ---------"<< cloudbox_limits[0]<<" "
               << cloudbox_limits[1] <<"\n";
          cout << endl;

          // Calculate abs_vec_array 
          abs_vec_array[p_index-cloudbox_limits[0]].resize(stokes_dim); 
          scat_p_index = p_index; 
          //(CE) p_index - cloudbox_limits[0] ????
          // This depends on the definition of pnd_field.
          abs_vec_agenda.execute();
          abs_vec_array[p_index-cloudbox_limits[0]] = abs_vec;
         
          // Calculate ext_mat_array
          ext_mat_array[p_index-cloudbox_limits[0]].
            resize(stokes_dim, stokes_dim); 
          scat_p_index = p_index;
          ext_mat_agenda.execute();
          ext_mat_array[p_index-cloudbox_limits[0]] = ext_mat;

           //Generate Planck function.
          Numeric T = t_field(p_index, 0, 0); 
          Numeric f = f_grid[scat_f_index];
          planck_vector[p_index-cloudbox_limits[0]] = planck(f, T);
         
          
        }//End of p_grid loop over the cloudbox
      

      //======================================================================
      // Radiative transfer inside the cloudbox
      //=====================================================================

      ext_mat.resize(stokes_dim, stokes_dim);
      abs_vec.resize(stokes_dim);
      sca_vec.resize(stokes_dim);
      stokes_vec.resize(stokes_dim);

      // Upward propagating direction, 90° is uplooking in the spherical 
      // geometry
      if(scat_za_grid[scat_za_index] <= 90.)
        {
          
          // Loop over all positions inside the cloud box defined by the 
          // cloudbox_limits exculding the upper boundary.
          for(Index p_index = cloudbox_limits[0]; p_index
                <= cloudbox_limits[1]-1; p_index ++)
            {
              
              //Initialize ppath for 1D.
              ppath_init_structure(ppath_step, 1, 1);
              // See documentation of ppath_init_structure for understanding
              // the parameters.
              
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
              
              // Length of the path between the two layers.
              l_step = ppath_step.l_step[0];
              
              // Check if the agenda has returned ppath.step with reasonable 
              // values. 
              cout << "\n ";
              PpathPrint( ppath_step, "ppath");
              
              
              // Calculate the average of the coefficients for the layers
              // to be considered in the 
              // radiative transfer calculation.
              
              for (Index i = 0; i < stokes_dim; i++)
                {
                  
                   // Extinction matrix requires a second loop over stokes_dim
                  for (Index j = 0; j < stokes_dim; j++)
                    {
                      ext_mat(i,j) = 0.5*( ext_mat_array
                                           [p_index-cloudbox_limits[0]](i,j) +
                                           ext_mat_array
                                           [p_index-cloudbox_limits[0]+1]
                                           (i,j));
                    }
                 
                  cout << "1.............."<< endl;
                  // Absorption vector
                  abs_vec[i] = 0.5*( abs_vec_array
                                     [p_index-cloudbox_limits[0]][i] +
                                     abs_vec_array
                                     [p_index-cloudbox_limits[0]+1][i] );
                 cout << "2.............."<< endl;
                  // Extract sca_vec from sca_field.
                  sca_vec[i] =
                    0.5*( scat_field(p_index-cloudbox_limits[0], 
                                     0, 0, scat_za_index, 0, i)
                          + scat_field(p_index-cloudbox_limits[0]+1,
                                       0, 0, scat_za_index, 0, i)) ;
                   cout << "3.............."<< endl;
                  // Extract stokes_vec from i_field.
                  stokes_vec[i] = 
                    0.5*( i_field((p_index-cloudbox_limits[0]), 0, 0,
                                  scat_za_index, 0, i)
                          +i_field((p_index-cloudbox_limits[0]) + 1,
				       0, 0, scat_za_index, 0, i));
                  
                }// Closes loop over stokes_dim.
              cout << "4.............."<< endl;
              //Planck function
              planck_function = 0.5*( planck_vector
                                     [p_index-cloudbox_limits[0]] +
                                     planck_vector
                                     [p_index-cloudbox_limits[0]+1]);

                   // Call scat_rte_agenda:
              scat_rte_agenda.execute();
              
              // Assign calculated Stokes Vector to i_field. 
              i_field((p_index - cloudbox_limits[0]), 0, 0, scat_za_index, 0,
                      Range(joker)) = stokes_vec;
              
            }// Close loop over p_grid (inside cloudbox).
        }// Closes 'if' for upward propagation direction.

          
      // Downward propagating direction:
      // There are some angles slightly above 90° where the same layer 
      // is intersected, this case has to be treated separately.
      else if(scat_za_grid[scat_za_index]  > 90.2)
        {
          // Loop over all positions inside the cloud box defined by the 
          // cloudbox_limits including the upper boundary.
          for(Index p_index = cloudbox_limits[0]+1; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
            
              //Initialize ppath for 1D.
              ppath_init_structure(ppath_step, 1, 1);
              // See documentation of ppath_init_structure for
              //understanding the
              // parameters.
              
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
              
              // Length of the path between the two layers.
              l_step = ppath_step.l_step[0];
              
              
              // Check if the agenda has returned ppath.step with reasonable 
              // values. 
              cout << "\n ";
              PpathPrint( ppath_step, "ppath");
              
              // Calculate averadge values of the coefficients.
 
              // Loop over the Stokes dimensions
              for (Index i = 0; i < stokes_dim; i++)
                {
                  // Extinction matrix requires a second loop over stokes_dim
                  for (Index j = 0; j < stokes_dim; j++)
                    {
                      ext_mat(i,j) = 0.5*( ext_mat_array
                                           [p_index-cloudbox_limits[0]](i,j) +
                                           ext_mat_array
                                           [p_index-cloudbox_limits[0]-1]
                                           (i,j));
                    }
                  // Absorption vector
                  abs_vec[i] = 0.5*( abs_vec_array
                                     [p_index-cloudbox_limits[0]][i] +
                                     abs_vec_array
                                     [p_index-cloudbox_limits[0]-1][i] );
                
                  // Extract sca_vec from sca_field.
                  sca_vec[i] =
                    0.5*( scat_field((p_index-cloudbox_limits[0]), 
                                     0, 0, scat_za_index, 0, i)
                          + scat_field((p_index-cloudbox_limits[0]) - 1,
                                       0, 0, scat_za_index, 0, i)) ;
                  // Extract stokes_vec from i_field.
                  stokes_vec[i] = 
                    0.5*( i_field((p_index-cloudbox_limits[0]), 0, 0,
                                  scat_za_index, 0, i)
                          +i_field((p_index-cloudbox_limits[0]) - 1,
				       0, 0, scat_za_index, 0, i));

                }// Closes loop over stokes_dim.  

              //Planck function
              planck_function = 0.5*( planck_vector
                                      [p_index-cloudbox_limits[0]] +
                                      planck_vector
                                      [p_index-cloudbox_limits[0]-1]);
                
              // Call scat_rte_agenda:
              scat_rte_agenda.execute();
              
              // Assign calculated Stokes Vector to i_field. 
              i_field((p_index - cloudbox_limits[0]), 0, 0, scat_za_index, 0,
                      Range(joker)) = stokes_vec;
            
            }// Close loop over p_grid (inside cloudbox).
        }// Closes 'if' for downward propagation direction.
      
      // FIXME: special case, close to 90°
      // now it takes the mean value
      // FIXMEEEEEEE !!!!!!!!!!!!!!!!!!!!
    //   else 
//         {
//           for(Index p_index = cloudbox_limits[0]; p_index
//                 <= cloudbox_limits[1]; p_index ++)
//             {
//               for (Index i = 0; i < stokes_dim; i++)
//                 { 
//                   Numeric i_field1 = i_field((p_index - cloudbox_limits[0]), 0,
//                                             0, scat_za_index-1, 0, i);
//                   Numeric i_field2 = i_field((p_index - cloudbox_limits[0]), 0,
//                                             0, scat_za_index+1, 0, i);
//                   i_field((p_index - cloudbox_limits[0]), 0, 0, scat_za_index,
//                           0, i) = 0.5 * (i_field1 + i_field2);
//                 } //end for stokes_dim
//             }// end for p_index
//         }//end special case
      
    } //Closes loop over scat_za_index.

} // End of the function.





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
  assert(stokes_dim <= 4 && stokes_dim > 0);
  
  stokes_vec.resize(stokes_dim);
 
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
  
\param scat_field Output : scattering integraal
\param pha_mat Output : phase matrix.
\param i_field Input : the intensity field 
\param pha_mat_spt Input : the phase matrix for single particle type
\param pnd_field Input : the particle number density field.
\param scat_za_grid zenith angle grid
\param scat_aa_grid azimuth angle grid
\param p_grid pressure grid
\param lat_grid latitude grid
\param lon_grid longitude grid
\param stokes_dim dimension of stokes vector
\param atmosphere_dim atmospheric dimension
\param cloudbox_limits Limits of the cloudbox.

\author Sreerekha Ravi
\date 2002-06-20
*/

void
scat_fieldCalc(//WS Output:
	       Tensor6& scat_field,
	       Tensor4& pha_mat,
	       //WS Input: 
	       const Tensor6& i_field,
	       const Tensor5& pha_mat_spt,
	       const Tensor4& pnd_field,
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
  Index N_pt = pha_mat_spt.nshelves();
  Index Nza = scat_za_grid.nelem();
  Index Naa = scat_aa_grid.nelem();
  Index Np  = i_field.nvitrines();//STR
  cout<<"Naa in the scattering integral routine"<<" "<<Naa<<"\n";
  Tensor4 product_field(Np,Nza, Naa, stokes_dim);//earlier tensor3; added pressure index STR
  Vector product_field_vec(stokes_dim);
  // scat_field is of the same size as *i_field*
  scat_field.resize( i_field.nvitrines(),
		     i_field.nshelves(),
		     i_field.nbooks(),
		     i_field.npages(),
		     i_field.nrows(),
		     i_field.ncols() );
  
  //Note that the size of i_field and scat_field are the same 
  if (atmosphere_dim ==3){
    
    assert ( is_size( i_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
		      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
		      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
		      scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(),
		      stokes_dim));
    assert ( is_size( scat_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
		      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
		      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
		      scat_za_grid.nelem(), 
		      scat_aa_grid.nelem(),
		      stokes_dim));
    
  }
  
  else if (atmosphere_dim == 1 ){
    assert ( is_size( i_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
		      1, 
		      1,
		      scat_za_grid.nelem(), 
		      1,
		      stokes_dim));
        assert ( is_size( scat_field, 
		      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
		      1, 
		      1,
		      scat_za_grid.nelem(), 
		      1,
		      stokes_dim));
  }
  //When atmospheric dimension , atmosphere_dim = 1
  if( atmosphere_dim == 1 ){
    
    /*there is a loop only over the pressure index when we calculate
      the pha_mat from pha_mat_spt and pnd_field using the method
      pha_matCalc.  Note that the
      lat_index and lon_index are set to zero*/
    
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++) 
      {
	pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
		    atmosphere_dim, p_index, 0, 
		    0);
      }
    
    
    // Get pha_mat at the grid positions
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++)
      {
	for (Index za_index = 0; za_index < Nza; ++ za_index)
	  {
	    for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
	      {
		
		ConstMatrixView pha = pha_mat(za_index,
					      aa_index,
					      Range(joker),
					      Range(joker));
		//cout<<"numer of vitrines"<<i_field.nvitrines()<<"\n";
		ConstVectorView ifield_in = 
		  i_field((p_index-cloudbox_limits[0]),
			  //changed p_index above
			  0,0,
			  za_index,
			  //aa_index,
			  0,
			  Range(joker));
		
		// just a multiplication of ifield_in and pha
		
		mult (product_field_vec, pha, ifield_in);
		
		for (Index i = 0; i < stokes_dim; i++)
		  {
		    
		    product_field((p_index-cloudbox_limits[0]),
				  //changed p_index
				  za_index,
				  aa_index, 
				  i) = product_field_vec[i];
		  }
	      }
	  }
      }
    
	    /*integration of the product of ifield_in and pha
	      over zenith angle and azimuth angle grid. It calls
	      here the integration routine AngIntegrate_trapezoid*/
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++)
      {
	for (Index za_index = 0; za_index < Nza; ++ za_index)
	  {
	    for (Index i = 0; i < stokes_dim; i++)
	      {
		MatrixView product_field_mat = product_field((p_index-cloudbox_limits[0]),//changed p_index
							     Range(joker),
							     Range(joker),
							     i);
		
	
		scat_field((p_index-cloudbox_limits[0]),//changed p_index
			   0,
			   0,
			   za_index, // by STR
			   //Range(joker),commented by STR
			   //Range(joker),commented by STR
			   0,
			   i)  =   AngIntegrate_trapezoid(product_field_mat,
							  scat_za_grid,
							  scat_aa_grid);
		
	      }
	  }//i closed it here 
	
      }
  }
  //cout<<"scat_field"<<scat_field<<"\n";	
  //exit(1);
  
  
  //When atmospheric dimension , atmosphere_dim = 3
  if( atmosphere_dim == 3 ){
    
    /*there is a loop only over the pressure index when we calculate
      the pha_mat from pha_mat_spt and pnd_field using the method
      pha_matCalc.  */
    
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++)
      {
	for (Index lat_index = cloudbox_limits[2]; lat_index <= 
	       cloudbox_limits[3]; lat_index++)
	  {
	    for (Index lon_index = cloudbox_limits[4]; lon_index <= 
		   cloudbox_limits[5]; lon_index++)
	      {
		pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
			    atmosphere_dim, p_index, lat_index, 
			    lon_index);
	      }
	  }
	
      }
    
    // Get pha_mat at the grid positions
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++)
      {
	for (Index lat_index = cloudbox_limits[2]; lat_index <= 
	       cloudbox_limits[3]; lat_index++)
	  {
	    for (Index lon_index = cloudbox_limits[4]; lon_index <= 
		   cloudbox_limits[5]; lon_index++)
	      {
		for (Index za_index = 0; za_index < Nza; ++ za_index)
		  {
		    for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
		      {			
			ConstMatrixView pha = pha_mat(za_index,
						      aa_index,
						      Range(joker),
						      Range(joker));
			
			//ifield_in is a vector of size 4 if stokes_dim =1
			ConstVectorView ifield_in =
			  i_field((p_index-cloudbox_limits[0]),
				  //changed p_index above
				  (lat_index-cloudbox_limits[2]),
				  //changed lat_index above
				  (lon_index-cloudbox_limits[4]),
				  //changed lat_index above
				  za_index,
				  aa_index,
				  Range(joker));
			
			// just a multiplication of ifield_in and pha
			
			mult (product_field_vec, pha, ifield_in);
			for (Index i = 0; i < stokes_dim; i++)
			  {
			    //??? should this inlcude latitude and longitude also 
			    product_field((p_index-cloudbox_limits[0]),
					  za_index,
					  aa_index,
					  i) = product_field_vec[i];
			  }
		      }
		  }
	      }
	  }
      }
    
    /*integration of the product of ifield_in and pha
      over zenith angle and azimuth angle grid. It 
      calls here the integration routine 
      AngIntegrate_trapezoid*/
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
	 p_index++)
      {
	for (Index lat_index = cloudbox_limits[2]; lat_index <= 
	       cloudbox_limits[3]; lat_index++)
	  {
	    for (Index lon_index = cloudbox_limits[4]; lon_index <= 
		   cloudbox_limits[5]; lon_index++)
	      {
		for (Index i = 0; i < stokes_dim; i++)
		  {
		    MatrixView product_field_mat = 
		      product_field((p_index-cloudbox_limits[0]),//changed p_index
				    Range(joker),
				    Range(joker),
				    i);
		    
		    scat_field((p_index-cloudbox_limits[0]),//changed p_index
			       lat_index,
			       lon_index,
			       Range(joker),
			       Range(joker),i)
		      = AngIntegrate_trapezoid(product_field_mat,
					       scat_za_grid,
					       scat_aa_grid);
		    
		  }
		
	      }
	  }
      }
  }
  
 
}
  

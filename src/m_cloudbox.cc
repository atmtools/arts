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
  \author Patrick Eriksson
  \date   2002-05-08 

  \brief  Workspace functions releated to the cloud box.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "arts.h"
#include "array.h"
#include "check_input.h"
#include "matpackI.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "cloudbox.h"
#include "logic.h"

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
  WS Input:
  \param ppath_step_agenda Agenda to compute a propagation path step.
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
  \param f_index       Frequency index.
  \param blackbody_ground Flag to treat ground as blackbody.
  \param stokes_dim    The number of Stokes components to be calculated.
  \param atmosphere_dim Atmospheric dimension.
*/
void i_fieldIterate(
		    // WS Output:
		    Tensor6& i_field,
		    Ppath& ppath_step,
		    Tensor6& i_field_old,
		    Tensor6& scat_field,
		    // WS Input:
		    const Agenda& ppath_step_agenda,
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

  // Copy i_field to i_field_old.
 
  i_field_old = i_field;
  
  //Calculate scattered field vector for all points in the cloudbox.

  scat_field = i_field;
  // ---- here will be the function to calculate the scattering integral
  
  
  //Update i_field.
  if( atmosphere_dim == 1 )
    {
      i_fieldUpdate1D(i_field, ppath_step, ppath_step_agenda, i_field_old, amp_mat, scat_field,
		       cloudbox_limits, scat_za_grid, scat_aa_grid, p_grid, 
		       lat_grid, lon_grid, t_field, z_field, z_ground,
		       r_geoid, f_grid, scat_f_index, blackbody_ground, 
		       stokes_dim);
    }
  
  //Convergence test has to be here.
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
  WS Input:
  \param ppath_step_agenda Agenda to compute a propagation path step.
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
*/
void i_fieldUpdate1D(// WS Output:
		     Tensor6& i_field,
		     Ppath& ppath_step,
		     // WS Input:
		     const Agenda& ppath_step_agenda,
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
	  Matrix ext_mat(4,4,3e-12);
  
	  //absorption vector
	  Vector abs_vec(4,1e-14);
	  
	  //scattering integral vector;
	  Vector sca_vec(4);
	  
	  //Extract sca_vec from sca_field.
	  for (Index i = 0; i < stokes_dim; i++)
	    {
	      sca_vec[i] = scat_field(p_index, 0, 0,
				      scat_za_index, 0, i);
	     }
	

	  //Generate Planck function.
	  Numeric T = t_field(p_index, 0, 0);
	  Numeric B;
	  Numeric f = f_grid[scat_f_index];
	  B = planck(f, T);

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
	  

	  // Call ppath_step_agenda 
	  // Check that the agenda takes te right input:

	  if( !ppath_step_agenda.is_output(ppath_step_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must generate ppath_step as an output.");
	    }

	  if( !ppath_step_agenda.is_input(atmosphere_dim_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use atmosphere_dim
                                   as an input.");
	    }

	  if( !ppath_step_agenda.is_input(p_grid_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use p_grid
                                   as an input.");
	    }
	  
	  if( !ppath_step_agenda.is_input(lat_grid_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use lat_grid
                                   as an input.");
	    }
	  
	  if( !ppath_step_agenda.is_input(lon_grid_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use lon_grid
                                   as an input.");
	    }

	  if( !ppath_step_agenda.is_input(z_field_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use z_field
                                   as an input.");
	    }
	  
	  if( !ppath_step_agenda.is_input(r_geoid_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use r_geoid
                                   as an input.");
	    }

	  if( !ppath_step_agenda.is_input(z_ground_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use z_ground
                                   as an input.");
	    }

	  if( !ppath_step_agenda.is_input(blackbody_ground_) )
	    {
	      throw runtime_error("The agenda ppath_step_agenda must use blackbody_ground
                                   as an input.");
	    }

	  
	  // Everything checked. Now call agenda.
	  ppath_step_agenda.execute();
	  
	  // Check if the agenda has returned ppath.step with reasonable values 
	  cout << "\n ";
	  PpathPrint( ppath_step, "ppath");
  

          Vector sto_vec(4);
	  // Perform RT calculation.
	  rte_scat_vecCalc(sto_vec, ext_mat, abs_vec, sca_vec, 
			   ppath_step.l_step[0], B); 
	  
	  // Assign claculated Stokes Vector to i_field. 
	  i_field(p_index, 0, 0, scat_za_index, 0, Range(joker)) = sto_vec;
  
	  // Close all loops.
	}
    }
}
		



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
  \author Patrick Eriksson, Claudia Emde and Sreerekha T. R.
  \date   2002-05-08 

  \brief  Workspace functions related to the definintion of the cloud box.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "xml_io.h"
#include "messages.h"
#include "gridded_fields.h"
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


//! This method interpolates clear sky field on the cloudbox boundary 
//on all grid points inside the cloud box. 

/*! 
  This method uses a linear 3D interpolation scheme to obtain the 
  radiation field on all grid points inside the cloud box form the 
  clear sky field on the cloud bod boundary.
  
  \param i_field Output : Intensity field
  \param scat_i_p Input : Intensity field on cloudbox boundary 
  (equal pressure surfaces)
  \param scat_i_lat Input : Intensity field on cloudbox boundary 
  (equal latitude surfaces)
  \param scat_i_lon Input : Intensity field on cloudbox boundary
  (equal longitude surfaces)
  \param f_grid Input : frequency grid
  \param scat_f_index Input : the frequency index for scattering calculation
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmospere_dim Input : dimension of atmosphere
*/
void i_fieldSetClearsky(Tensor6& i_field,
		const Tensor7& scat_i_p,
		const Tensor7& scat_i_lat,
		const Tensor7& scat_i_lon,
		const Vector& f_grid,
		const Index& scat_f_index,
		const Vector& p_grid,
		const Vector& lat_grid,
		const Vector& lon_grid,
		const ArrayOfIndex& cloudbox_limits,
		const Index& atmosphere_dim )
{
 
  if(atmosphere_dim == 1)
    {
      Index  N_f = scat_i_p.nlibraries();
      if (f_grid.nelem() != N_f){
	
	throw runtime_error(" scat_i_p should have same frequency  "
			    " dimension as f_grid");
      }
     
      if(scat_i_p.nvitrines() != 2){
	throw runtime_error("scat_i_p should have only two elements "
			    "in pressure grid which corresponds "
			    "to the two pressure surfaces");
      }
      
     
      Index N_za = scat_i_p.npages() ;
   
      Index N_aa = scat_i_p.nrows();
   
      Index N_i = scat_i_p.ncols();
         
      //1. interpolation - pressure grid
      
      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1,
		     1,
		     1,
		     N_za,
		     N_aa,
		     N_i);
      
      /*the old grid is having only two elements, corresponding to the 
	cloudbox_limits and the new grid have elements corresponding to
	all grid points inside the cloudbox plus the cloud_box_limits*/

      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      
      gridpos(p_gp,
	      p_grid[Range(cloudbox_limits[0], 
			   2,
			   (cloudbox_limits[1]- cloudbox_limits[0]))],
	      p_grid[Range(cloudbox_limits[0], 
			   (cloudbox_limits[1]- cloudbox_limits[0])+1)]);
      
      Matrix itw((cloudbox_limits[1]- cloudbox_limits[0])+1, 2);
      interpweights ( itw, p_gp );
      
      for (Index za_index = 0; za_index < N_za ; ++ za_index)
	{
	  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
	    {
	      for (Index i = 0 ; i < N_i ; ++ i)
		{
		  
		  VectorView target_field = i_field(Range(joker),
						    0,
						    0,
						    za_index,
						    aa_index,
						    i);
		  
		  ConstVectorView source_field = scat_i_p(scat_f_index,
							  Range(joker),    
							  0,
							  0,
							  za_index,
							  aa_index,
							  i);
		  
		  interp(target_field,
			 itw,
			 source_field,
			 p_gp);
		}
	      
	    }
	}

    }
  if(atmosphere_dim == 3)
    {
      Index  N_f = scat_i_p.nlibraries();
      if (scat_i_lat.nlibraries() != N_f || 
	  scat_i_lon.nlibraries() != N_f){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same frequency dimension");
      }
      Index N_p = scat_i_lat.nvitrines();
      if(scat_i_lon.nvitrines() != N_p ||
	 p_grid.nelem()         != N_p ){
	throw runtime_error("scat_i_lat and scat_i_lon should have  "
			    "same pressure grid dimension as p_grid");
      }
      
      Index N_lat = scat_i_p.nshelves();
      
      if(scat_i_lon.nshelves() != N_lat ||
	 lat_grid.nelem()      != N_lat){
	throw runtime_error("scat_i_p and scat_i_lon should have  "
			    "same latitude grid dimension as lat_grid");
      }
  
      Index N_lon = scat_i_p.nbooks();
      if(scat_i_lat.nbooks() != N_lon ||
	 lon_grid.nelem()    != N_lon ){
	throw runtime_error("scat_i_p and scat_i_lat should have  "
			    "same longitude grid dimension as lon_grid");
      }
      if(scat_i_p.nvitrines() != 2){
	throw runtime_error("scat_i_p should have only two elements "
			    "in pressure grid which corresponds "
			    "to the two pressure surfaces");
      }
      
      if(scat_i_lat.nshelves() != 2){
	throw runtime_error("scat_i_lat should have only two elements "
			    "in latitude grid which corresponds "
			    "to the two latitude surfaces");
	
      }
      if(scat_i_lon.nbooks() != 2){
	throw runtime_error("scat_i_lon should have only two elements "
			    "in longitude grid which corresponds "
			    "to the two longitude surfaces");
	
      }
      Index N_za = scat_i_p.npages() ;
      if (scat_i_lat.npages() != N_za || 
	  scat_i_lon.npages() != N_za){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same dimension for zenith angles");
      }
      Index N_aa = scat_i_p.nrows();
      if (scat_i_lat.nrows() != N_aa || 
	  scat_i_lon.nrows() != N_aa){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same dimension for azimuth angles");
      }
      Index N_i = scat_i_p.ncols();
      if (scat_i_lat.ncols() != N_i || 
	  scat_i_lon.ncols() != N_i){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same value for stokes_dim and can take only"
			    "values 1,2,3 or 4");
      }
      
      //1. interpolation - pressure grid, latitude grid and longitude grid
      
   
      //i_field
      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
		     (cloudbox_limits[3]- cloudbox_limits[2])+1,
		     (cloudbox_limits[5]- cloudbox_limits[4])+1,
		     N_za, 
		     N_aa,
		     N_i);
      
      
      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      ArrayOfGridPos lat_gp((cloudbox_limits[3]- cloudbox_limits[2])+1);
      ArrayOfGridPos lon_gp((cloudbox_limits[5]- cloudbox_limits[4])+1);

      /*the old grid is having only two elements, corresponding to the 
	cloudbox_limits and the new grid have elements corresponding to
	all grid points inside the cloudbox plus the cloud_box_limits*/
      
      gridpos(p_gp,
	      p_grid[Range(cloudbox_limits[0], 
			   2,
			   (cloudbox_limits[1]- cloudbox_limits[0]))],
	      p_grid[Range(cloudbox_limits[0], 
			   (cloudbox_limits[1]- cloudbox_limits[0])+1)]);
      gridpos(lat_gp,
	      lat_grid[Range(cloudbox_limits[2], 
			     2,
			     (cloudbox_limits[3]- cloudbox_limits[2]))],
	      lat_grid[Range(cloudbox_limits[2], 
			     (cloudbox_limits[3]- cloudbox_limits[2])+1)]);
      gridpos(lon_gp,
	      lon_grid[Range(cloudbox_limits[4], 
			     2,
			     (cloudbox_limits[5]- cloudbox_limits[4]))],
	      lon_grid[Range(cloudbox_limits[4], 
			     (cloudbox_limits[5]- cloudbox_limits[4])+1)]);
      
      //interpolation weights corresponding to pressure, latitude and 
      //longitude grids.

      Matrix itw_p((cloudbox_limits[1]- cloudbox_limits[0])+1, 2);
      Matrix itw_lat((cloudbox_limits[3]- cloudbox_limits[2])+1, 2);
      Matrix itw_lon((cloudbox_limits[5]- cloudbox_limits[4])+1, 2);

      interpweights ( itw_p, p_gp );
      interpweights ( itw_lat, lat_gp );
      interpweights ( itw_lon, lon_gp );

      // interpolation - pressure grid
      for (Index lat_index = cloudbox_limits[2]; 
	   lat_index < cloudbox_limits[3] ; ++ lat_index)
	{
	  for (Index lon_index = cloudbox_limits[4]; 
	       lon_index < cloudbox_limits[5] ; ++ lon_index)
	    {
	      for (Index za_index = 0; za_index < N_za ; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
		    {
		      for (Index i = 0 ; i < N_i ; ++ i)
			{
			  
			  VectorView target_field = i_field(Range(joker),
							    lat_index,
							    lon_index,
							    za_index,
							    aa_index,
							    i);
			  
			  ConstVectorView source_field = scat_i_p(scat_f_index,
								  Range(joker),    
								  lat_index,
								  lon_index,
								  za_index,
								  aa_index,
								  i);
			  
			  interp(target_field,
				 itw_p,
				 source_field,
				 p_gp);
			}
		    }
		}
	    } 
	}
      //interpolation latitude
      for (Index p_index = cloudbox_limits[0]; 
	   p_index < cloudbox_limits[1] ; ++ p_index)
	{
	  for (Index lon_index = cloudbox_limits[4]; 
	       lon_index < cloudbox_limits[5] ; ++ lon_index)
	    {
	      for (Index za_index = 0; za_index < N_za ; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
		    {
		      for (Index i = 0 ; i < N_i ; ++ i)
			{
			  
			  VectorView target_field = i_field(p_index,
							    Range(joker),
							    lon_index,
							    za_index,
							    aa_index,
							    i);
			  
			  ConstVectorView source_field = scat_i_p(scat_f_index,
								  p_index,
								  Range(joker),    
								  lon_index,
								  za_index,
								  aa_index,
								  i);
			  
			  interp(target_field,
				 itw_lat,
				 source_field,
				 lat_gp);
			}
		    }
		}
	    } 
	}
      //interpolation -longitude
      for (Index p_index = cloudbox_limits[0]; 
	   p_index < cloudbox_limits[1] ; ++ p_index)
	{
	  for (Index lat_index = cloudbox_limits[2]; 
	       lat_index < cloudbox_limits[3] ; ++ lat_index)
	    {
	      for (Index za_index = 0; za_index < N_za ; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
		    {
		      for (Index i = 0 ; i < N_i ; ++ i)
			{
			  
			  VectorView target_field = i_field(p_index,
							    lat_index,
							    Range(joker),
							    za_index,
							    aa_index,
							    i);
			  
			  ConstVectorView source_field = scat_i_p(scat_f_index,
								  p_index,    
								  lat_index,
								  Range(joker),
								  za_index,
								  aa_index,
								  i);
			  
			  interp(target_field,
				 itw_lon,
				 source_field,
				 lon_gp);
			}
		    }
		}
	    } 
	}
      //end of interpolation
    }//ends atmosphere_dim = 3
}


      
//! Set constant initial field

/*! 
  This method sets the initial field inside the cloudbox to a constant value.

  The user can specify a value for each Stokes dimension in the control file
  by the variable i_field_value, which is a vector containing 4 elements, the
  value of the initial field for each Stokes dimension.
  
  \param i_field Output : Intensity field
  \param scat_i_p Input : Intensity field on cloudbox boundary 
  (equal pressure surfaces)
  \param scat_i_lat Input : Intensity field on cloudbox boundary 
  (equal latitude surfaces)
  \param scat_i_lon Input : Intensity field on cloudbox boundary
  (equal longitude surfaces)
  \param f_grid Input : frequency grid
  \param scat_f_index Input : the frequency index for scattering calculation
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmospere_dim Input : dimension of atmosphere
  \param i_field_value Keyword: value of the constant initial field 

  \author Claudia Emde
  \date 2002-08-26
  
*/
void i_fieldSetConst(//WS Output:
                        Tensor6& i_field,
                        //WS Input:
                        const Tensor7& scat_i_p,
                        const Tensor7& scat_i_lat,
                        const Tensor7& scat_i_lon,
                        const Vector& f_grid,
                        const Index& scat_f_index,
                        const Vector& p_grid,
                        const Vector& lat_grid,
                        const Vector& lon_grid,
                        const ArrayOfIndex& cloudbox_limits,
                        const Index& atmosphere_dim,
                        const Index& stokes_dim,
                        // Keyword       
                        const Vector& i_field_values)
{
  if(atmosphere_dim == 1)
  {
    // In the 1D case the atmospheric layers are defined by p_grid and the
    // required interface is scat_i_p.
    Index N_za = scat_i_p.npages();
    Index N_p = p_grid.nelem();
   
    // Define the size of i_field.
    i_field.resize((cloudbox_limits[1] - cloudbox_limits[0])+1, 1, 1,  N_za,
                   1, 1);

    // Loop over all zenith angle directions.
    for (Index za_index = 0; za_index < N_za; za_index++)
      {
        for (Index i = 0; i < stokes_dim; i++)
          { 
            //set the value for the upper boundary
            i_field(cloudbox_limits[1]-cloudbox_limits[0], 0, 0, za_index,
                    0, i) = 
              scat_i_p(0, 1, 0, 0, za_index, 0, i);
            //set the value for the lower boundary 
            i_field(0, 0, 0, za_index, 0, i) = 
              scat_i_p(0, 0, 0, 0, za_index, 0, i);
            for (Index scat_p_index = 1; scat_p_index < cloudbox_limits[1] - 
                   cloudbox_limits[0]; scat_p_index++ )
              // The field inside the cloudbox is set to some arbitrary value.
              i_field(scat_p_index, 0, 0, za_index, 0, i) =  i_field_values[i];
          }    
      }
  }
  
  if(atmosphere_dim == 3)
    {
      throw runtime_error(
                          "i_fieldSetConst for 3D atmosphere will be \n"
                          "implemented soon."
                          );
    }
}



//! Initialize variables containing information about the particles.
/*! 
 
  WS Output:
  \param amp_mat_raw    Amplitude matrix data.
  \param pnd_field_raw  Particle number density field data.

*/
void ParticleTypeInit( //WS Output:
                      ArrayOfArrayOfTensor6& amp_mat_raw,
                      ArrayOfArrayOfTensor3& pnd_field_raw
                      )
{
  amp_mat_raw.reserve(20);
  pnd_field_raw.reserve(20); 
}


//! Read amplitute matrix and particle number density field from data base.
/*! 
  This method allows the user to chose particle types and particle number 
  density fields. 
  There is one database for particle number density fields ( ....),
  which includes the following particle types:

  Another database (....) containes the amplitude matrices for those particle
  types from which all optical properties can be derived.
  
  \param amp_mat_raw Amplitude matrix data.
  \param pnd_field_raw Particle number density field data.
  \param amp_mat_file Filename for amplitude matrix data.
  \param pnd_field_file Filename for pnd field data.
*/
void ParticleTypeAdd( //WS Output:
                 ArrayOfArrayOfTensor6& amp_mat_raw,
                 ArrayOfArrayOfTensor3& pnd_field_raw,
                 // Keyword:
                 const String& amp_mat_file,
                 const String& pnd_field_file)
{

  // Append *amp_mat_raw* and *pnd_field_raw* with empty Arrays of Tensors. 
  ArrayOfTensor6 amp_mat_data;
  amp_mat_raw.push_back(amp_mat_data);
  
  //ArrayOfTensor3 pnd_field_data;
  //pnd_field_raw.push_back(pnd_field_data);
    
  // Read amplitude matrix data:
  out2 << "Read amplitude matrix_data\n";
  read_gridded_tensor6( amp_mat_file, amp_mat_raw[amp_mat_raw.nelem()-1]);
  //read_gridded_tensor3 (pnd_field_file, pnd_field_raw[pnd_field_raw.nelem()-1]);
       
}


//! Method for the communication between cloudbox and clearsky.
/*
  This method puts the scattered radiation field into the interface variables 
  between the cloudbox and the clearsky, which are *scat_i_p*, *scat_i_lat*,
  *scat_i_lon*. As i_field is only stored for one frequency given by
  *scat_f_index* this method has to be executed after each scattering 
  calculation to store the scattered field on the boundary of the cloudbox.
  
  The best way to calculate spectra including the influence of 
  scattering is to set up the *scat_mono_agenda* where this method 
  can be included.
 
  
 \param scat_i_p i_field on pressure boundaries
 \param scat_i_lat i_field on latitude boundaries
 \param scat_i_lon i_field on longitude boundaries
 \param i_field radiation field inside cloudbox
 \param f_grid frequency grid
 \param scat_f_index index for scattering calculation
 \param p_grid pressure grid
 \param lat_grid latitude grid
 \param lon_grid longitude grid
 \param scat_za_grid zenith angle grid
 \param scat_aa_grid azimuth angle grid
 \param stokes_dim stokes dimension
 \param atmosphere_dim atmospheric dimension
 \param cloudbox_limits limits of the cloudbox

 \author Claudia Emde
 \date 2002-09-09
 
*/
void scat_iPut(//WS Output:
               Tensor7&  scat_i_p,
               Tensor7& scat_i_lat,
               Tensor7& scat_i_lon,
               //WS Input:
               const Tensor6& i_field,
               const Vector& f_grid,
               const Index& scat_f_index,
               const Vector& p_grid,
               const Vector& lat_grid,
               const Vector& lon_grid,
               const Vector& scat_za_grid,
               const Vector& scat_aa_grid,
               const Index& stokes_dim,
               const Index& atmosphere_dim,
               const ArrayOfIndex& cloudbox_limits)
{
  // Some sizes:
  Index N_f = f_grid.nelem();
  Index N_p = p_grid.nelem();
  Index N_lat = lat_grid.nelem();
  Index N_lon = lon_grid.nelem();
  Index N_za = scat_za_grid.nelem();
  Index N_aa = scat_aa_grid.nelem();

  // Some checks:
  assert ( is_size( i_field, 
                    (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                    1, 
                    1,
                    N_za, 
                    1,
                    stokes_dim));

  assert ( is_size( scat_i_p,
                    N_f, 2, N_lat, N_lon, N_za, N_aa, stokes_dim ));

  assert ( is_size( scat_i_lat,
                    N_f, N_p, 2, N_lon, N_za, N_aa, stokes_dim ));
  
  assert ( is_size( scat_i_lon,
                    N_f, N_p, N_lat, 2, N_za, N_aa, stokes_dim ));

  // Put the i_field at the cloudbox boundary into the interface variable 
  // scat_i_p.
  if(atmosphere_dim == 1)
    {
      for (Index za = 0; za < N_za; za++)
            {
              for (Index i = 0; i < stokes_dim; i++)
                {  
                  //i_field at lower boundary
                  scat_i_p(scat_f_index, 0, 0, 0, za, 0, i) = 
                    i_field(0, 0, 0, za, 0, i);
                  //i_field at upper boundary
                  scat_i_p(scat_f_index, 1, 0, 0, za, 0, i) = 
                    i_field(cloudbox_limits[1] - cloudbox_limits[0],
                            0, 0, za, 0, i); 
                  // For 1D scat_i_lat and scat_i_lon are 0.
                  for(Index j = 0; j<2; j++)
                    {
                      scat_i_lat(scat_f_index, j, 0, 0, za, 0, i) = 0.;
                      scat_i_lon(scat_f_index, j, 0, 0, za, 0, i) = 0.;
                    }
                }//end stokes_dim
            }//end za loop
    }//end atmosphere_dim = 1
                                                        
      

      
  if(atmosphere_dim == 3)
    {
      throw runtime_error(
                          "i_fieldSetConst for 3D atmosphere will be"
                          "implemented soon."
                          );
    }
}


//! Scattered radiance on the cloudbox boundary.
/* 
 This method returns the radiances for a given direction and position on the 
 boundary of the cloudbox. It interpolates from *scat_za_grid* on the 
 requested direction. The variable *y_scat* is a matrix with the 
 dimensions [f_grid, stokes_dim].
  
  \param y_scat Scattered radiance.
  \param scat_i_p i_field on pressure boundaries.
  \param scat_i_lat i_field on latitude boundaries.
  \param scat_i_lon i_field on longitude boundaries.
  \param cloudbox_pos Position on the cloudbox boundary.
  \param cloudbox_los Direction of radiation.
  \param cloudbox_limits Cloudbox limits.
  \param atmosphere_dim Atmospheric dimension.
  \param stokes_dim Stokes dimension.

  \author Claudia Emde
  \date 2002-09-10

 */    
void y_scatCalc(//WS Output:
                Matrix& y_scat,
                //WS Input:
                const Tensor7& scat_i_p,
                const Tensor7& scat_i_lat,
                const Tensor7& scat_i_lon,
                const Vector& cloudbox_pos,
                const Vector& cloudbox_los,
                const ArrayOfIndex& cloudbox_limits,
                const Index& atmosphere_dim,
                const Index& stokes_dim,
                const Vector& scat_za_grid,
                const Vector& scat_aa_grid,
                const Vector& f_grid)
{

 if(atmosphere_dim == 1)
   {
     if (cloudbox_pos[0] != cloudbox_limits[0] &&
         cloudbox_pos[0] != cloudbox_limits[1])
       throw runtime_error(
                           "*cloudbox_pos* has to be on the boundary of the "
                           "cloudbox defined by *cloudbox_limits*."
                           );
     
     //Define a vector to interpolate the outgoing radiance which is 
     //defined on scat_za_grid on the requested zenith angle in 
     //*cloudbox_los*.
     Vector zenith_angle(1);
     zenith_angle[0] = cloudbox_los[0];
         
     //Array to store grid positions
     ArrayOfGridPos gp(1);
     //Matrix to store interpolation weights
     Matrix itw(scat_za_grid.nelem(),2);

     for(Index i = 0; i < stokes_dim; i++)
       {
         for(Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
           {
             //This vvariable holds the radiation for a specified frequency.
             //It is neccessairy because the interpolation is done for 
             //each frequency separately.
             Vector y_scat_f(scat_za_grid.nelem());

             //lower boundary
             if(cloudbox_pos[0] == cloudbox_limits[0])
               {
                 ConstVectorView y_f = scat_i_p(f_index, 0, 0, 0, 
                                                 Range(joker), 0, i);
                 y_scat_f = y_f;
               }
             //upper boundary
             else if(cloudbox_pos[0] == cloudbox_limits[1])
               {
                 ConstVectorView y_f = scat_i_p(f_index, 1, 0, 0,
                                                 Range(joker), 0, i);
                 y_scat_f = y_f;
               }
             //Define vector for the interpolated radiance.
             Vector y_scat_los(1);
             
             //Do the interpolation:
             interp(y_scat_los, itw, y_scat_f, gp);
             
             //Put the value into the matrix:
             y_scat(f_index, i) = y_scat_los[0];
           }//end frequency loop
       }//end stokes_dim loop
   }// end atmosphere_dim 1


 if(atmosphere_dim == 3)
    {
      throw runtime_error(
                          "i_fieldSetConst for 3D atmosphere will be \n"
                          "implemented soon."
                          );
    }
   
}

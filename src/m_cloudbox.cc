/* Copyright (C) 2002,2003
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
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
#include <cstdlib>
#include <cmath>

#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "xml_io.h"
#include "messages.h"
#include "gridded_fields.h"
#include "logic.h"
#include "rte.h"
#include "interpolation.h"
#include "special_interp.h"
#include "cloudbox.h"
#include "optproperties.h"
#include "math_funcs.h"
#include "physics_funcs.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! cloudboxOff
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void cloudboxOff(
        // WS Output:
        Index&           cloudbox_on,
        ArrayOfIndex&    cloudbox_limits,
        Agenda&          iy_cloudbox_agenda )
{
  cloudbox_on = 0;
  cloudbox_limits.resize(0);
  iy_cloudbox_agenda.resize(0);
}


//! cloudboxSetEmpty
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Claudia Emde
   \date   2002-05-11
*/
void cloudboxSetEmpty(//WS Output:
                      Tensor4& pnd_field,
                      ArrayOfSingleScatteringData& scat_data_raw,
                      //WS Input:
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid)
{
  // 3D  atmosphere
  if (lat_grid.nelem()>0)
    {
      //Resize pnd_field and set it to 0:
      pnd_field.resize(1, p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
      pnd_field = 0.;
    }
  else // 1D atmosphere
     {
      //Resize pnd_field and set it to 0:
      pnd_field.resize(1, p_grid.nelem(), 1, 1);
      pnd_field = 0.;
     }
  
  //Resize scat_data_raw and set it to 0:
  // Number iof particle types
  scat_data_raw.resize(1);
  scat_data_raw[0].ptype = PTYPE_MACROS_ISO;
  scat_data_raw[0].description = " ";
  // Grids which contain full ranges which one wants to calculate
  nlinspace(scat_data_raw[0].f_grid, 1e9, 10000e9, 5);  
  nlinspace(scat_data_raw[0].T_grid, 0, 400, 5);
  nlinspace(scat_data_raw[0].za_grid, 0, 180, 5);
  nlinspace(scat_data_raw[0].aa_grid, 0, 360, 5);
  // Resize the data arrays
  scat_data_raw[0].pha_mat_data.resize(5,5,5,1,1,1,6);
  scat_data_raw[0].pha_mat_data = 0.;
  scat_data_raw[0].ext_mat_data.resize(5,5,1,1,1);
  scat_data_raw[0].ext_mat_data = 0.;
  scat_data_raw[0].abs_vec_data.resize(5,5,1,1,1);
  scat_data_raw[0].abs_vec_data = 0.;
}



//! cloudboxSetManually
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-19
*/
void cloudboxSetManually(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
        Index&          scat_za_interp,
        // WS Input:
        const Index&    atmosphere_dim,
        const Vector&   p_grid,
        const Vector&   lat_grid,
        const Vector&   lon_grid,
        // Control Parameters:
        const Numeric& p1,
        const Numeric& p2,
        const Numeric& lat1,
        const Numeric& lat2,
        const Numeric& lon1,
        const Numeric& lon2 )
{
  // Default interpolation is linear (scat_za_interp = 0)
  scat_za_interp = 0;

  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  if( atmosphere_dim == 2 )
    { throw runtime_error( "The cloud box is not defined for 2D." ); }

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
      for( cloudbox_limits[0]=1; p_grid[cloudbox_limits[0]+1]>=p1; 
                                                     cloudbox_limits[0]++ ) {}
    }
  if( p2 < p_grid[p_grid.nelem()-2] )
    {
      cloudbox_limits[1] = p_grid.nelem() - 1;
    }
  else
    {
      for( cloudbox_limits[1]=p_grid.nelem()-2; 
                    p_grid[cloudbox_limits[1]-1]<=p2; cloudbox_limits[1]-- ) {}
    }

  // Latitude limits
  if( atmosphere_dim >= 2 )
    {
      for( cloudbox_limits[2]=1; lat_grid[cloudbox_limits[2]+1]<=lat1; 
                                                     cloudbox_limits[2]++ ) {}
      for( cloudbox_limits[3]=lat_grid.nelem()-2; 
                lat_grid[cloudbox_limits[3]-1]>=lat2; cloudbox_limits[3]-- ) {}
    }

  // Longitude limits
  if( atmosphere_dim == 3 )
    {
      for( cloudbox_limits[4]=1; lon_grid[cloudbox_limits[4]+1]<=lon1; 
                                                     cloudbox_limits[4]++ ) {}
      for( cloudbox_limits[5]=lon_grid.nelem()-2; 
                lon_grid[cloudbox_limits[5]-1]>=lon2; cloudbox_limits[5]-- ) {}
    }
}



//! cloudboxSetManuallyAltitude
/*!
  This is basically the same function as *cloudboxSetManually* but one 
  can specify the cloudbox altitude instead of pressure.
  It takes pressure indices corresponding to the requested cloud altitude 
  at the first latitude/longitude limit.

  See the the online help (arts -d FUNCTION_NAME)

  \author Claudia Emde
  \date   2004-02-20
*/
void cloudboxSetManuallyAltitude(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
        Index&          scat_za_interp,
        // WS Input:
        const Index&    atmosphere_dim,
        const Tensor3&  z_field,
        const Vector&   lat_grid,
        const Vector&   lon_grid,
        // Control Parameters:
        const Numeric& z1,
        const Numeric& z2,
        const Numeric& lat1,
        const Numeric& lat2,
        const Numeric& lon1,
        const Numeric& lon2 )
{
  // Default interpolation is linear (scat_za_interp = 0)
  scat_za_interp = 0;

  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  if( atmosphere_dim == 2 )
    { throw runtime_error( "The cloud box is not defined for 2D." ); }

  // Check keyword arguments
  if( z1 >= z2 )
    throw runtime_error( 
                        "The altitude in *z1* must be smaller than the altitude in *z2*." );
  if( z1 <= z_field(0, 0, 0) )
    throw runtime_error( "The altitude in *z1* must be larger than the "
                         "first value in *z_field*." );
  if( z2 >= z_field(z_field.npages()-1, 0, 0) )
    throw runtime_error( "The altitude in *z2* must be smaller than the "
                         "last value in *z_field*." );
  if( atmosphere_dim == 3 )
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
      if( lon2 <= lon1 )
        throw runtime_error
          ("The longitude in *lon2* must be bigger than the longitude in *lon1*.");
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

  // Pressure/altitude limits
  if( z1 < z_field(1, 0, 0) )
    {
      cloudbox_limits[0] = 0;
    }
  else
    {
      for( cloudbox_limits[0]=1; z_field(cloudbox_limits[0]+1, 0, 0) <= z1; 
                                                     cloudbox_limits[0]++ ) {}
    }
  if( z2 > z_field(z_field.npages()-2, 0, 0) )
    {
      cloudbox_limits[1] = z_field.npages() - 1;
    }
  else
    {
      for( cloudbox_limits[1]=z_field.npages()- 2; 
           z_field(cloudbox_limits[1]-1, 0, 0) >= z2; cloudbox_limits[1]-- ) {}
    }

  // Latitude limits
  if( atmosphere_dim >= 2 )
    {
      for( cloudbox_limits[2]=1; lat_grid[cloudbox_limits[2]+1]<=lat1; 
                                                     cloudbox_limits[2]++ ) {}
      for( cloudbox_limits[3]=lat_grid.nelem()-2; 
                lat_grid[cloudbox_limits[3]-1]>=lat2; cloudbox_limits[3]-- ) {}
    }

  // Longitude limits
  if( atmosphere_dim == 3 )
    {
      for( cloudbox_limits[4]=1; lon_grid[cloudbox_limits[4]+1]<=lon1; 
                                                     cloudbox_limits[4]++ ) {}
      for( cloudbox_limits[5]=lon_grid.nelem()-2; 
                lon_grid[cloudbox_limits[5]-1]>=lon2; cloudbox_limits[5]-- ) {}
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
  \param f_index Input : the frequency index for scattering calculation
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmosphere_dim Input : dimension of atmosphere
  
  \author Sreerekha T.R.
   \date   2002-07-24
*/
void i_fieldSetClearsky(Tensor6& i_field,
                const Tensor7& scat_i_p,
                const Tensor7& scat_i_lat,
                const Tensor7& scat_i_lon,
                const Vector& f_grid,
                const Index& f_index,
                const Vector& p_grid,
                const Vector& lat_grid,
                const Vector& lon_grid,
                const ArrayOfIndex& cloudbox_limits,
                const Index& atmosphere_dim )
{
  
  out2 << "Interpolate boundary clearsky field to obtain the initial field.\n";

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

      i_field = 0.;

      

      /*the old grid is having only two elements, corresponding to the 
        cloudbox_limits and the new grid have elements corresponding to
        all grid points inside the cloudbox plus the cloud_box_limits*/

      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      
      p2gridpos(p_gp,
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
                  
                  ConstVectorView source_field = scat_i_p(f_index,
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
      Index N_p = cloudbox_limits[1] - cloudbox_limits[0] + 1;
      if(scat_i_lon.nvitrines() != N_p ||
         scat_i_lat.nvitrines() != N_p ){
        throw runtime_error("scat_i_lat and scat_i_lon should have  "
                            "same pressure grid dimension as p_grid");
      }
      
      Index N_lat =  cloudbox_limits[3] - cloudbox_limits[2] + 1;
      
      if(scat_i_lon.nshelves() != N_lat ||
         scat_i_p.nshelves()   != N_lat){
        throw runtime_error("scat_i_p and scat_i_lon should have  "
                            "same latitude grid dimension as lat_grid");
      }
  
      Index N_lon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      if(scat_i_lat.nbooks() != N_lon ||
         scat_i_p.nbooks()   != N_lon ){
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
      
      p2gridpos(p_gp,
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
      for (Index lat_index = 0; 
           lat_index <= (cloudbox_limits[3]-cloudbox_limits[2]); ++ lat_index)
        {
          for (Index lon_index = 0; 
               lon_index <= (cloudbox_limits[5]-cloudbox_limits[4]);
                            ++ lon_index)
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
                          
                          ConstVectorView source_field = scat_i_p(f_index,
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
      for (Index p_index = 0; 
           p_index <= (cloudbox_limits[1]-cloudbox_limits[0]) ; ++ p_index)
        {
          for (Index lon_index = 0; 
               lon_index <= (cloudbox_limits[5]-cloudbox_limits[4]) ;
               ++ lon_index)
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
                          
                          ConstVectorView source_field = scat_i_lat(f_index,
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
      for (Index p_index = 0; 
           p_index <= (cloudbox_limits[1]-cloudbox_limits[0]); ++ p_index)
        {
          for (Index lat_index = 0; 
               lat_index <= (cloudbox_limits[3]-cloudbox_limits[2]);
               ++ lat_index)
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
                          
                          ConstVectorView source_field = scat_i_lon(f_index,
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
  The method works only for monochromatic calculations (number of elements 
  in f_grid =1).

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
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmosphere_dim Input : dimension of atmosphere
  \param stokes_dim     FIXME: Add documentation.
  \param i_field_values Keyword : value of the constant initial field 

  \author Claudia Emde
  \date 2002-08-26
  
*/
void i_fieldSetConst(//WS Output:
                        Tensor6& i_field,
                        //WS Input:
                        const Tensor7& scat_i_p,
                        const Tensor7& scat_i_lat,
                        const Tensor7& scat_i_lon,
                        const Vector& p_grid,
                        const Vector& lat_grid,
                        const Vector& lon_grid,
                        const ArrayOfIndex& cloudbox_limits,
                        const Index& atmosphere_dim,
                        const Index& stokes_dim,
                        // Keyword       
                        const Vector& i_field_values)
{
  out2 << "Set initial field to constant values: " << i_field_values << "\n"; 

  // In the 1D case the atmospheric layers are defined by p_grid and the
  // required interface is scat_i_p.
  Index N_za = scat_i_p.npages();
  Index N_aa = scat_i_p.nrows();
  Index N_i = stokes_dim; 
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Grids have to be adapted to atmosphere_dim.
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  
   // Check the input:
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*"); 


 
  if(atmosphere_dim == 1)
    {
      cout << "atm_dim = 1" << endl; 
      
    // Define the size of i_field.
    i_field.resize((cloudbox_limits[1] - cloudbox_limits[0])+1, 1, 1,  N_za,
                   1, N_i);
    i_field = 0.;

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
  else 
    {
      if ( !is_size(scat_i_p, 1, 2, Nlat_cloud, 
                Nlon_cloud, N_za, N_aa, stokes_dim)  
           || !is_size(scat_i_lat, 1, Np_cloud, 2, 
                       Nlon_cloud, N_za, N_aa, stokes_dim)  
           || !is_size(scat_i_lon, 1, Np_cloud,  
                       Nlat_cloud, 2, N_za, N_aa, stokes_dim) )
        throw runtime_error(
                            "One of the interface variables (*scat_i_p*, "
                            "*scat_i_lat* or *scat_i_lon*) does not have "
                            "the right dimensions.  \n Probably you have "
                            "calculated them before for another value of "
                            "*stokes_dim*."
                            );
      

      
      cout << "atm_dim = 3" << endl;      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
                     (cloudbox_limits[3]- cloudbox_limits[2])+1,
                     (cloudbox_limits[5]- cloudbox_limits[4])+1,
                     N_za, 
                     N_aa,
                     N_i);
      
      i_field = 0.;
      

      // Loop over all directions:
      for (Index za_index = 0; za_index < N_za; za_index++)
        {
          for (Index aa_index = 0; aa_index < N_aa; aa_index++)
            {
              for (Index i = 0; i < stokes_dim; i++)
                { 
                  // pressure boundaries
                  // 
                   for (Index lat_index = cloudbox_limits[2]; 
                         lat_index <= cloudbox_limits[3]; lat_index++)
                      {
                        for (Index lon_index = cloudbox_limits[4]; 
                             lon_index <= cloudbox_limits[5]; lon_index++)
                          {
                            //set the value for the upper pressure boundary
                            i_field(cloudbox_limits[1]-cloudbox_limits[0], 
                                    lat_index-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_p(0, 1, lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i);
                            //set the value for the lower pressure boundary 
                            i_field(0, lat_index-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) =  
                              scat_i_p(0, 0, lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i);
                          }
                      }
                   
                   for (Index p_index = cloudbox_limits[0]; 
                        p_index <= cloudbox_limits[1]; p_index++)
                     {
                       // latitude boundaries
                       //
                        for (Index lon_index = cloudbox_limits[4]; 
                             lon_index <= cloudbox_limits[5]; lon_index++)
                          {
                            // first boundary
                            i_field(p_index-cloudbox_limits[0], 
                                    cloudbox_limits[3]-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_lat(0, p_index-cloudbox_limits[0],
                                         1, lon_index-cloudbox_limits[4],
                                         za_index, aa_index, i);
                            // second boundary
                            i_field(p_index-cloudbox_limits[0], 0, 
                                    lon_index-cloudbox_limits[4], 
                                    za_index, aa_index, i) =  
                              scat_i_lat(0, p_index-cloudbox_limits[0], 0,
                                         lon_index-cloudbox_limits[4],
                                         za_index, aa_index, i);
                            
                          }                                                    
                        // longitude boundaries
                   for (Index lat_index = cloudbox_limits[2]; 
                             lat_index <= cloudbox_limits[3]; lat_index++)
                          {
                            // first boundary
                            i_field(p_index-cloudbox_limits[0],
                                    lat_index-cloudbox_limits[2],
                                    cloudbox_limits[5]-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_lon(0, p_index-cloudbox_limits[0],
                                         lat_index-cloudbox_limits[2], 1,
                                         za_index, aa_index, i);
                            // second boundary
                            i_field(p_index-cloudbox_limits[0],  
                                    lat_index-cloudbox_limits[2],
                                    0, 
                                    za_index, aa_index, i) =  
                              scat_i_lon(0, p_index-cloudbox_limits[0],
                                         lat_index-cloudbox_limits[2], 0,
                                         za_index, aa_index, i);
                          } //lat_grid loop
                     } //p_grid loop
                   //
                   // Set the initial field to a constant value inside the 
                   // cloudbox:
                   // 
                   for( Index p_index = (cloudbox_limits[0]+1); 
                                 p_index <  cloudbox_limits[1] ;
                                 p_index ++)
                     {
                       for (Index lat_index = (cloudbox_limits[2]+1); 
                            lat_index < cloudbox_limits[3]; 
                            lat_index++)
                         {
                           for (Index lon_index = (cloudbox_limits[4]+1); 
                                lon_index < cloudbox_limits[5];
                                lon_index++)
                             {
                               i_field(p_index-cloudbox_limits[0],
                                       lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i) =  
                                 i_field_values[i];
                             }
                         }
                     }
                } // stokes loop
            } // aa_grid loop
        } // za_grid loop
       
    } // atmosphere dim = 3
}


//! Initialize variables containing information about the particles.
/*! 
 
  WS Output:
  \param scat_data_raw  Single scattering data.
  \param pnd_field_raw  Particle number density field data.

*/
void ParticleTypeInit( //WS Output:
                      ArrayOfSingleScatteringData& scat_data_raw,
                      ArrayOfGriddedField3& pnd_field_raw
                      )
{
  scat_data_raw.reserve(20);
  pnd_field_raw.reserve(20); 
}


//! Read single scattering data and particle number densities
/*! 
 The particle number densities *pnd_field_raw* can also be generated 
 without using the ARTS database. In this case it is easier to generate one 
 file including the data for all particle types than generating 
 a pnd file for each particle type. This means, that for such cases it is 
 more handy to us *ParticleTypeAddAll* instead of *ParticleTypeAdd*.
  
  \param scat_data_raw Single scattering data.
  \param pnd_field_raw Particle number density field data.
  \param filename_scat_data Filename for a file containing the filenames of the 
                       single scattering data files.
  \param pnd_field_file Filename for pnd field data.

  \author Claudia Emde
  \date 2004-03-08

  */
void ParticleTypeAddAll( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfGriddedField3&  pnd_field_raw,
                 // Keyword:
                 const String& filename_scat_data,
                 const String& pnd_field_file)
{
 
  
  ArrayOfString data_files;
  xml_read_from_file(filename_scat_data, data_files);

  scat_data_raw.resize(data_files.nelem());
  
  for (Index i = 0; i<data_files.nelem(); i++)
    {
      
      out2 << "Read single scattering data\n";
      xml_read_from_file( data_files[i], 
                          scat_data_raw[i]);
    }
  
  out2 << "Read particle number density date \n";
  xml_read_from_file(pnd_field_file, pnd_field_raw);
       
}


//! Read single scattering data and particle number density field from 
//  data base
/*! 
  This method allows the user to chose hydro-meteor species and particle
  number density fields. The data is added to *scat_data_raw* and 
  *pnd_field_raw*. 
  
  There is one database for particle number density fields ( ....),
  which includes the following particle types:

  Another database (....) contains the single scattering properties for 
  hydro-meteor species.
  
  \param scat_data_raw Single scattering data.
  \param pnd_field_raw Particle number density field data.
  \param scat_data_file Filename for scattering data.
  \param pnd_field_file Filename for pnd field data.

  \author Claudia Emde
  \date 2003-02-24

  \date 2004-02-23 Included structure ArrayOfGriddedField3.
*/
void ParticleTypeAdd( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfGriddedField3&  pnd_field_raw,
                 // Keyword:
                 const String& scat_data_file,
                 const String& pnd_field_file)
{
  
  // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors. 
  SingleScatteringData scat_data;
  scat_data_raw.push_back(scat_data);
  
  GriddedField3 pnd_field_data;
  pnd_field_raw.push_back(pnd_field_data);
  
  out2 << "Read single scattering data\n";
  
  xml_read_from_file( scat_data_file, scat_data_raw[scat_data_raw.nelem()-1]);
  
  out2 << "Read particle number density date \n";
  if (pnd_field_file.nelem()>0)
    {
     xml_read_from_file(pnd_field_file,pnd_field_raw[pnd_field_raw.nelem()-1]);
    }   
   
}


//! Calculate particle number density fields.
/*!
  This method interpolates the data for particle number density fields
  on the atmospheric grids used for the calculation.

   See also the the online help (arts -d FUNCTION_NAME)

   \author Sreerekha T.R.
   \date   2003-04-17

   \date   2004-02-23 (CE:) Used ArrayOfGriddedField3 instead of 
           ArrayOfArrayOfTensor3. 

*/

void pnd_fieldCalc(//WS Output:
		   Tensor4& pnd_field,
                   //WS Input
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
		   const ArrayOfGriddedField3& pnd_field_raw,
                   const Index& atmosphere_dim
                   )
{

  // Basic checks of input variables
  //
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  
  //==========================================================================
  if ( atmosphere_dim == 1)
    {
      
      //Resize variables
      pnd_field.resize(pnd_field_raw.nelem(), p_grid.nelem(), 1, 1);
      
      // Gridpositions:
      ArrayOfGridPos gp_p(p_grid.nelem());
         
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, pnd_field_raw[i].p_grid, p_grid);
      
          // Interpolation weights:
	   Matrix itw(p_grid.nelem(), 2);
	   // (2 interpolation weights are required for 1D interpolation)
          interpweights( itw, gp_p);
               // Interpolate:
          interp( pnd_field(i, joker, 0, 0),
                  itw, pnd_field_raw[i].data(joker, 0, 0), gp_p);
        }
      
    }

  //=========================================================================
  else if(atmosphere_dim == 2)
    {
      //Resize variables
      pnd_field.resize(pnd_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       1);
      
      
      // Gridpositions:
      ArrayOfGridPos gp_p(p_grid.nelem());
      ArrayOfGridPos gp_lat(lat_grid.nelem());
                  
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, pnd_field_raw[i].p_grid, p_grid);
          gridpos(gp_lat, pnd_field_raw[i].data(0, joker, 0), 
                  lat_grid);
	  
          // Interpolation weights:
	  Tensor3 itw(p_grid.nelem(), lat_grid.nelem(), 4);
	  // (8 interpolation weights are required for 3D interpolation)
          interpweights( itw, gp_p, gp_lat);
          
          // Interpolate:
          interp( pnd_field(i, joker, joker, 0),
                  itw, pnd_field_raw[i].data(joker, joker, 0),
                  gp_p, gp_lat);
        }
    }

  //================================================================
  // atmosphere_dim = 3    
  else
    {
      //Resize variables
      pnd_field.resize(pnd_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       lon_grid.nelem());
      
      
      // Gridpositions:
      ArrayOfGridPos gp_p(p_grid.nelem());
      ArrayOfGridPos gp_lat(lat_grid.nelem());
      ArrayOfGridPos gp_lon(lon_grid.nelem());
      
      
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, pnd_field_raw[i].p_grid, p_grid);
          gridpos(gp_lat, pnd_field_raw[i].lat_grid, lat_grid);
          gridpos(gp_lon, pnd_field_raw[i].lon_grid, lon_grid);
          
          // Interpolation weights:
          Tensor4 itw(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem(), 8);
          // (8 interpolation weights are required for 3D interpolation)
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( pnd_field(i, joker, joker, joker),
                  itw, pnd_field_raw[i].data, gp_p, gp_lat, gp_lon);
        }
    }
}


//! Method for the communication between cloudbox and clearsky.
/*
  This method puts the scattered radiation field into the interface variables 
  between the cloudbox and the clearsky, which are *scat_i_p*, *scat_i_lat*,
  *scat_i_lon*. As i_field is only stored for one frequency given by
  *f_index* this method has to be executed after each scattering 
  calculation to store the scattered field on the boundary of the cloudbox.
  
  The best way to calculate spectra including the influence of 
  scattering is to set up the *scat_mono_agenda* where this method 
  can be included.
 
  
 \param scat_i_p i_field on pressure boundaries
 \param scat_i_lat i_field on latitude boundaries
 \param scat_i_lon i_field on longitude boundaries
 \param i_field radiation field inside cloudbox
 \param f_grid frequency grid
 \param f_index index for scattering calculation
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
               const Index& f_index,
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
  Index N_za = scat_za_grid.nelem();
  Index N_aa = scat_aa_grid.nelem();
  Index N_p = cloudbox_limits[1] - cloudbox_limits[0] +1;

  // Some checks:
  
  assert( f_index < f_grid.nelem() );
 
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Grids have to be adapted to atmosphere_dim.
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  
   // Check the input:
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*"); 
  // End of checks.

  // Put the i_field at the cloudbox boundary into the interface variable 
  // scat_i_p.
  if(atmosphere_dim == 1)
    {
      // Check size of i_field.
      assert ( is_size( i_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                        1, 
                        1,
                        N_za, 
                        1,
                        stokes_dim));
      
      assert ( is_size( scat_i_p,
                        N_f, 2, 1, 1, N_za, 1, stokes_dim ));
      
      for (Index za = 0; za < N_za; za++)
        {
          for (Index i = 0; i < stokes_dim; i++)
            {  
              
              //i_field at lower boundary
              scat_i_p(f_index, 0, 0, 0,
                       za, 0, i) = 
                i_field(0, 0, 0, za, 0, i);
              //i_field at upper boundary
              scat_i_p(f_index, 1, 0, 0,
                       za, 0, i) = 
                i_field(cloudbox_limits[1] - cloudbox_limits[0],
                        0, 0, za, 0, i); 
              
            }//end stokes_dim
        }//end za loop
    }//end atmosphere_dim = 1
        
  if(atmosphere_dim == 3)
    {
      // Some sizes relevant for 3D atmosphere
      Index N_lat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      Index N_lon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      
      // Check size of i_field.
      assert ( is_size( i_field, 
                        cloudbox_limits[1] - cloudbox_limits[0] + 1,
                        N_lat,
                        N_lon,
                        N_za, 
                        N_aa,
                        stokes_dim));

      // Resize interface variables:
      scat_i_p.resize(N_f, 2, N_lat, N_lon, N_za, N_aa, stokes_dim);
      scat_i_lat.resize(N_f, N_p, 2, N_lon, N_za, N_aa, stokes_dim);
      scat_i_lon.resize(N_f, N_p, N_lat, 2, N_za, N_aa, stokes_dim);
 
      for (Index za = 0; za < N_za; za++)
        {
          for (Index aa = 0; aa < N_aa; aa++)
            {
              for (Index i = 0; i < stokes_dim; i++)
                {  
                  //
                  // Put i_field in scat_i_p:
                  //
                  for (Index lat = 0; lat < N_lat; lat++)
                    {
                      for (Index lon = 0; lon < N_lon; lon++)
                        {
                          //i_field at lower boundary
                          scat_i_p(f_index, 0, lat, lon,
                                   za, aa, i) = 
                            i_field(0, lat, lon, za, aa, i);
                          //i_field at upper boundary
                          scat_i_p(f_index, 1, lat, lon,
                                   za, aa, i) = 
                            i_field(cloudbox_limits[1]-cloudbox_limits[0],
                                    lat, lon, za, aa, i);
                        }
                    }
                  // 
                  // Put i_field in scat_i_lat:
                  //
                  for (Index p = 0; p < N_p; p++)
                    {
                      for (Index lon = 0; lon < N_lon; lon++)
                        {
                          //i_field at lower boundary
                          scat_i_lat(f_index, p, 0, lon,
                                     za, aa, i) = 
                            i_field(p, 0, lon, za, aa, i);
                          //i_field at upper boundary
                          scat_i_lat(f_index, p, 1, lon,
                                     za, aa, i) = 
                            i_field(p, cloudbox_limits[3]-
                                    cloudbox_limits[2],
                                    lon, za, aa, i);
                        }
                      //
                      // Put i_field in scat_i_lon:
                      for (Index lat = 0; lat < N_lat; lat++)
                        {
                          //i_field at lower boundary
                          scat_i_lon(f_index, p, lat, 0,
                                     za, aa, i) = 
                            i_field(p, lat, 0, za, aa, i);
                          //i_field at upper boundary
                          scat_i_lon(f_index, p, lat, 1,
                                     za, aa, i) = 
                            i_field(p, lat, cloudbox_limits[5]-
                                    cloudbox_limits[4], za, aa, i);
                        } 
                    }
                }
            }
        }
    }
}



//! CloudboxGetIncoming
/*! 
   See the the online help (arts -d FUNCTION_NAME)

  \author Sreerekha T.R., Claudia Emde
  \date 2002-10-07
  \date 2004-09-27   Modified and adapted to iy_calc by Patrick Eriksson.
*/
void CloudboxGetIncoming(
              Tensor7&        scat_i_p,
              Tensor7&        scat_i_lat,
              Tensor7&        scat_i_lon,
              Matrix&         iy,
              Ppath&          ppath,
              Ppath&          ppath_step,
              Vector&         rte_pos,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Vector&         rte_los,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         iy_surface_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Vector&         scat_za_grid,
        const Vector&         scat_aa_grid )
{
  Index Nf       = f_grid.nelem();
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index Nza      = scat_za_grid.nelem();
  Index Ni       = stokes_dim;


  //--- Check input ----------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  if( cloudbox_on == 0  ||  cloudbox_limits.nelem() == 0 )
    throw runtime_error( "The cloudbox must be activated, and it is not.");
  if( scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180. )
        throw runtime_error(
                     "*scat_za_grid* must include 0° and 180° as endpoints." );
  //--------------------------------------------------------------------------


  // Dummy variable for flag cloudbox_on. It has to be 0 here not to get
  // stuck in an infinite loop (if some propagation path hits the cloud
  // box at some other position.
  Index cloudbox_on_dummy = 0;

  // Make all agendas silent
  const Index   agenda_verb = true;


  if( atmosphere_dim == 1 )
    {
      // Resize interface variables:
      scat_i_p.resize( Nf, 2, 1, 1, Nza, 1, Ni );
      scat_i_lat.resize( 0, 0, 0, 0, 0, 0, 0 );
      scat_i_lon.resize( 0, 0, 0, 0, 0, 0, 0 );

      //Define the variables for position and direction.
      Vector   los(1), pos(1);

      //--- Get scat_i_p at lower boundary
      pos[0] = r_geoid(0,0) + z_field( cloudbox_limits[0], 0, 0 );

      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index ++)
        {
          los[0] =  scat_za_grid[scat_za_index];

          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat,
                   rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
                   iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
                   atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                   r_geoid, z_surface, cloudbox_on_dummy,  cloudbox_limits, 
                   pos, los, f_grid, stokes_dim, agenda_verb );

          scat_i_p( joker, 0, 0, 0, scat_za_index, 0, joker ) = iy;
        }

      //--- Get scat_i_p at upper boundary
      pos[0] = r_geoid(0,0) + z_field( cloudbox_limits[1], 0, 0 );

      for (Index scat_za_index = 0; scat_za_index < Nza;  scat_za_index ++)
        {
          los[0] =  scat_za_grid[scat_za_index];

          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat,
                   rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
                   iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
                   atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                   r_geoid, z_surface, cloudbox_on_dummy,  cloudbox_limits, 
                   pos, los, f_grid, stokes_dim, agenda_verb );

          scat_i_p( joker, 1, 0, 0, scat_za_index, 0, joker ) = iy;
        }
    }


  //--- atmosphere_dim = 3: --------------------------------------------------
  else
    {
      Index Naa = scat_aa_grid.nelem();

      if( scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360. )
        throw runtime_error(
                     "*scat_aa_grid* must include 0° and 360° as endpoints." );

      Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      
       // Convert scat_za_grid to "sensor coordinates"
      //(-180° < azimuth angle < 180°)
      //
      Vector aa_grid(Naa);
      for(Index i = 0; i<Naa; i++)
        aa_grid[i] = scat_aa_grid[i] - 180;

      // Resize interface variables:
      scat_i_p.resize( Nf, 2, Nlat_cloud, Nlon_cloud, Nza, Naa, Ni );
      scat_i_lat.resize( Nf, Np_cloud, 2, Nlon_cloud, Nza, Naa, Ni );
      scat_i_lon.resize( Nf, Np_cloud, Nlat_cloud, 2, Nza, Naa, Ni );

      // Define the variables for position and direction.
      Vector   los(2), pos(3);

      
      //--- Get scat_i_p at lower boundary
      for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                               lon_index + cloudbox_limits[4]) 
                     + z_field(cloudbox_limits[0],
                               lat_index + cloudbox_limits[2],
                               lon_index + cloudbox_limits[4]);
              pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
              pos[2] = lon_grid[lon_index + cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      los[0] = scat_za_grid[scat_za_index];
                      los[1] = aa_grid[scat_aa_index];
                      
                      // For end points of scat_za_index, we need only to
                      // perform calculations for first scat_aa
                      if( !( ( scat_za_index == 0  ||  
                               scat_za_index == (Nza-1) )  &&  
                             scat_aa_index == 0 ) )
                        {
                          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, 
                          rte_gp_lat, rte_gp_lon, rte_los, ppath_step_agenda, 
                          rte_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cloudbox_agenda,  atmosphere_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                          cloudbox_on_dummy,  cloudbox_limits, 
                          pos, los, f_grid, stokes_dim, agenda_verb );
                        }

                      scat_i_p( joker, 0, lat_index, lon_index, 
                                scat_za_index, scat_aa_index, joker) = iy;
                    }
                }
            }
        }
      

      //--- Get scat_i_p at upper boundary
      for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                               lon_index + cloudbox_limits[4]) 
                     + z_field(cloudbox_limits[1],
                               lat_index + cloudbox_limits[2],
                               lon_index + cloudbox_limits[4]);
              pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
              pos[2] = lon_grid[lon_index + cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      los[0] = scat_za_grid[scat_za_index];
                      los[1] = aa_grid[scat_aa_index];
                      
                      // For end points of scat_za_index, we need only to
                      // perform calculations for first scat_aa
                      if( !( ( scat_za_index == 0  ||  
                               scat_za_index == (Nza-1) )  &&  
                             scat_aa_index == 0 ) )
                        {
                          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, 
                          rte_gp_lat, rte_gp_lon, rte_los, ppath_step_agenda, 
                          rte_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cloudbox_agenda,  atmosphere_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                          cloudbox_on_dummy,  cloudbox_limits, 
                          pos, los, f_grid, stokes_dim, agenda_verb );
                        }

                      scat_i_p( joker, 1, lat_index, lon_index, 
                                scat_za_index, scat_aa_index, joker) = iy;
                    }
                }
            }
        }
              
       
      //--- Get scat_i_lat (1st boundary):
      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              pos[0] = r_geoid(cloudbox_limits[2],
                               lon_index + cloudbox_limits[4]) 
                     + z_field(p_index + cloudbox_limits[0],
                               cloudbox_limits[2], 
                               lon_index + cloudbox_limits[4]);
              pos[1] = lat_grid[cloudbox_limits[2]];
              pos[2] = lon_grid[lon_index + cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      los[0] = scat_za_grid[scat_za_index];
                      los[1] = aa_grid[scat_aa_index];
                      
                      // For end points of scat_za_index, we need only to
                      // perform calculations for first scat_aa
                      if( !( ( scat_za_index == 0  ||  
                               scat_za_index == (Nza-1) )  &&  
                             scat_aa_index == 0 ) )
                        {
                          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, 
                          rte_gp_lat, rte_gp_lon, rte_los, ppath_step_agenda, 
                          rte_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cloudbox_agenda,  atmosphere_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                          cloudbox_on_dummy,  cloudbox_limits, 
                          pos, los, f_grid, stokes_dim, agenda_verb );
                        }

                      scat_i_lat( joker, p_index, 0, lon_index,
                                  scat_za_index, scat_aa_index, joker) = iy;
                    }
                }
            }
        }

      
      //--- Get scat_i_lat (2nd boundary)
      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              pos[0] = r_geoid(cloudbox_limits[3],
                               lon_index + cloudbox_limits[4])
                     + z_field(p_index + cloudbox_limits[0],
                               cloudbox_limits[3],
                               lon_index + cloudbox_limits[4]);
              pos[1] = lat_grid[cloudbox_limits[3]];
              pos[2] = lon_grid[lon_index + cloudbox_limits[4]];
               
              for (Index scat_za_index = 0; scat_za_index < Nza;
                    scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      los[0] = scat_za_grid[scat_za_index];
                      los[1] = aa_grid[scat_aa_index];
                      
                      // For end points of scat_za_index, we need only to
                      // perform calculations for first scat_aa
                      if( !( ( scat_za_index == 0  ||  
                               scat_za_index == (Nza-1) )  &&  
                             scat_aa_index == 0 ) )
                        {
                          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, 
                          rte_gp_lat, rte_gp_lon, rte_los, ppath_step_agenda, 
                          rte_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cloudbox_agenda,  atmosphere_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                          cloudbox_on_dummy,  cloudbox_limits, 
                          pos, los, f_grid, stokes_dim, agenda_verb );
                        }

                      scat_i_lat( joker, p_index, 1, lon_index, 
                                  scat_za_index, scat_aa_index, joker) = iy;
                    }
                }
            }
        }    


      //--- Get scat_i_lon (1st boundary):
      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
            {
              pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                               cloudbox_limits[4]) 
                     + z_field(p_index + cloudbox_limits[0],
                               lat_index + cloudbox_limits[2],
                               cloudbox_limits[4]);
              pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
              pos[2] = lon_grid[cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      los[0] = scat_za_grid[scat_za_index];
                      los[1] = aa_grid[scat_aa_index];
                      
                      // For end points of scat_za_index, we need only to
                      // perform calculations for first scat_aa
                      if( !( ( scat_za_index == 0  ||  
                               scat_za_index == (Nza-1) )  &&  
                             scat_aa_index == 0 ) )
                        {
                          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, 
                          rte_gp_lat, rte_gp_lon, rte_los, ppath_step_agenda, 
                          rte_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cloudbox_agenda,  atmosphere_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                          cloudbox_on_dummy,  cloudbox_limits, 
                          pos, los, f_grid, stokes_dim, agenda_verb );
                        }

                      scat_i_lon( joker, p_index, lat_index, 0, 
                                  scat_za_index, scat_aa_index, joker) = iy;
                    }
                }
            }
        }

      
      //--- Get scat_i_lon (2nd boundary)
      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
            {
              pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                               cloudbox_limits[5]) 
                     + z_field(p_index + cloudbox_limits[0],
                               lat_index + cloudbox_limits[2],
                               cloudbox_limits[5]);
              pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
              pos[2] = lon_grid[cloudbox_limits[5]];
              
              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      los[0] = scat_za_grid[scat_za_index];
                      los[1] = aa_grid[scat_aa_index];
                      
                      // For end points of scat_za_index, we need only to
                      // perform calculations for first scat_aa
                      if( !( ( scat_za_index == 0  ||  
                               scat_za_index == (Nza-1) )  &&  
                             scat_aa_index == 0 ) )
                        {
                          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, 
                          rte_gp_lat, rte_gp_lon, rte_los, ppath_step_agenda, 
                          rte_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cloudbox_agenda,  atmosphere_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                          cloudbox_on_dummy,  cloudbox_limits, 
                          pos, los, f_grid, stokes_dim, agenda_verb );
                        }

                      scat_i_lon( joker, p_index, lat_index, 1, 
                                  scat_za_index, scat_aa_index, joker) = iy;
                    }
                }
            }
        }
    }// End atmosphere_dim = 3.
}



//! CloudboxGetIncoming1DAtm
/*! 
   See the the online help (arts -d FUNCTION_NAME)

  \author Sreerekha T.R., Claudia Emde
  \date 2002-10-07
  \date 2004-10-01   Modified and adapted to iy_calc by Patrick Eriksson.
*/
void CloudboxGetIncoming1DAtm(
              Tensor7&        scat_i_p,
              Tensor7&        scat_i_lat,
              Tensor7&        scat_i_lon,
              Matrix&         iy,
              Ppath&          ppath,
              Ppath&          ppath_step,
              Vector&         rte_pos,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Vector&         rte_los,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         iy_surface_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Vector&         scat_za_grid,
        const Vector&         scat_aa_grid )
{
  Index Nf       = f_grid.nelem();
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
  Index Nza      = scat_za_grid.nelem();
  Index Naa      = scat_aa_grid.nelem();
  Index Ni       = stokes_dim;


  //--- Check input ----------------------------------------------------------
  if( atmosphere_dim != 3 )
    throw runtime_error( "The atmospheric dimensionality must be 3.");
  if( cloudbox_on == 0  ||  cloudbox_limits.nelem() == 0 )
    throw runtime_error( "The cloudbox must be activated, and it is not.");
  if( scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180. )
    throw runtime_error(
                     "*scat_za_grid* must include 0° and 180° as endpoints." );
  if( scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360. )
    throw runtime_error(
                     "*scat_aa_grid* must include 0° and 360° as endpoints." );
  //--------------------------------------------------------------------------

  // Dummy variable for flag cloudbox_on. It has to be 0 here not to get
  // stuck in an infinite loop (if some propagation path hits the cloud
  // box at some other position.
  Index cloudbox_on_dummy = 0;

  // Make all agendas silent
  const Index   agenda_verb = true;

      
  // Convert scat_za_grid to "sensor coordinates"
  //(-180° < azimuth angle < 180°)
  //
  Vector aa_grid(Naa);
  for(Index i = 0; i<Naa; i++)
    aa_grid[i] = scat_aa_grid[i] - 180;

  // As the atmosphere is spherically symmetric we only have to calculate 
  // one azimuth angle.
  Index aa_index = 0;

  // Resize interface variables:
  scat_i_p.resize(Nf, 2, Nlat_cloud, Nlon_cloud, Nza, Naa, Ni);
  scat_i_lat.resize(Nf, Np_cloud, 2, Nlon_cloud, Nza, Naa, Ni);
  scat_i_lon.resize(Nf, Np_cloud, Nlat_cloud, 2, Nza, Naa, Ni);
  
  // Define the variables for position and direction.
  Vector   los(2), pos(3);

  Index lat_index = 0;
  Index lon_index = 0;

  // These variables are constant for all calculations:
  pos[1] = lat_grid[cloudbox_limits[2]];
  pos[2] = lon_grid[lon_index + cloudbox_limits[4]];
  los[1] = aa_grid[aa_index];

  // Get scat_i_p at lower boundary
  pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                   lon_index + cloudbox_limits[4]) 
          + z_field(cloudbox_limits[0],
                    lat_index + cloudbox_limits[2],
                    lon_index + cloudbox_limits[4]);

  for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
    {
      los[0] = scat_za_grid[scat_za_index];
      
      iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
               rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
               iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
               z_surface, cloudbox_on_dummy,  cloudbox_limits, 
               pos, los, f_grid, stokes_dim, agenda_verb );

      for (Index lat = 0; lat < Nlat_cloud; lat ++)
        {
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_p( joker, 0, lat, lon, scat_za_index, aa, joker )
                    = iy;
                }
            }
        }
    }
      
  // Get scat_i_p at upper boundary
  pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                   lon_index + cloudbox_limits[4]) 
         + z_field(cloudbox_limits[1],
                   lat_index + cloudbox_limits[2],
                   lon_index + cloudbox_limits[4]);
  
  for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
    {
      los[0] = scat_za_grid[scat_za_index];
      
      iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
               rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
               iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
               z_surface, cloudbox_on_dummy,  cloudbox_limits, 
               pos, los, f_grid, stokes_dim, agenda_verb );

      for (Index lat = 0; lat < Nlat_cloud; lat ++)
        {
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_p( joker, 1, lat, lon, scat_za_index, aa, joker ) 
                    = iy;
                }
            }
        }
    }

              
       
  // Get scat_i_lat (1st boundary):
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      pos[0] = r_geoid(cloudbox_limits[2],
                       lon_index + cloudbox_limits[4]) 
             + z_field(p_index + cloudbox_limits[0],
                       cloudbox_limits[2],
                       lon_index + cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
        {
          los[0] = scat_za_grid[scat_za_index];
          
          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
               rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
               iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
               z_surface, cloudbox_on_dummy,  cloudbox_limits, 
               pos, los, f_grid, stokes_dim, agenda_verb );

          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lat( joker, p_index, 0, lon, scat_za_index, aa, joker)
                    = iy;
                }
            }
        }
    }
    
  // Get scat_i_lat (2nd boundary)
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      pos[0] = r_geoid(cloudbox_limits[3],
                       lon_index + cloudbox_limits[4]) 
             + z_field(p_index + cloudbox_limits[0],
                       cloudbox_limits[3],
                       lon_index + cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
        {
          los[0] = scat_za_grid[scat_za_index];
          
          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
               rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
               iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
               z_surface, cloudbox_on_dummy,  cloudbox_limits, 
               pos, los, f_grid, stokes_dim, agenda_verb );
          
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lat( joker, p_index, 1, lon, scat_za_index, aa, joker)
                    = iy;
                }
            }
        }
    }    

  // Get scat_i_lon (1st boundary):
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                       cloudbox_limits[4]) 
             + z_field(p_index + cloudbox_limits[0],
                       lat_index + cloudbox_limits[2],
                       cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
        {
          los[0] = scat_za_grid[scat_za_index];
          
          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
               rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
               iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
               z_surface, cloudbox_on_dummy,  cloudbox_limits, 
               pos, los, f_grid, stokes_dim, agenda_verb );

          for (Index lat = 0; lat < Nlat_cloud; lat ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lon( joker, p_index, lat, 0, scat_za_index, aa, joker)
                    = iy;
                }
            }
        }
    }
  
  // Get scat_i_lon (2nd boundary)
    for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                       cloudbox_limits[5]) 
             + z_field(p_index + cloudbox_limits[0],
                       lat_index + cloudbox_limits[2],
                       cloudbox_limits[5]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
        {
          los[0] = scat_za_grid[scat_za_index];
          
          iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
               rte_gp_lon, rte_los, ppath_step_agenda, rte_agenda, 
               iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda, 
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
               z_surface, cloudbox_on_dummy,  cloudbox_limits, 
               pos, los, f_grid, stokes_dim, agenda_verb );

          for (Index lat = 0; lat < Nlat_cloud; lat ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lon( joker, p_index, lat, 1, scat_za_index, aa, joker)
                    = iy;
                }
            }
        }
    }
}



//! iyInterpCloudboxField
/*! 
   See the the online help (arts -d FUNCTION_NAME)

  \author Patrick Eriksson
  \date 2004-09-29
*/
void iyInterpCloudboxField(
            Matrix&         iy,
      const Tensor7&        scat_i_p,
      const Tensor7&        scat_i_lat,
      const Tensor7&        scat_i_lon,
      const GridPos&        rte_gp_p,
      const GridPos&        rte_gp_lat,
      const GridPos&        rte_gp_lon,
      const Vector&         rte_los,
      const Index&          cloudbox_on,
      const ArrayOfIndex&   cloudbox_limits,
      const Index&          atmosphere_dim,
      const Index&          stokes_dim,
      const Vector&         scat_za_grid,
      const Vector&         scat_aa_grid,
      const Vector&         f_grid )
{
  iy_interp_cloudbox_field( iy, scat_i_p, scat_i_lat, scat_i_lon, rte_gp_p, 
                            rte_gp_lat, rte_gp_lon, rte_los, cloudbox_on, 
                            cloudbox_limits, atmosphere_dim, stokes_dim, 
                            scat_za_grid, scat_aa_grid, f_grid, "linear" );
}



//! iyInterpCubicCloudboxField
/*! 
   See the the online help (arts -d FUNCTION_NAME)

  \author Patrick Eriksson
  \date 2004-09-29
*/
void iyInterpCubicCloudboxField(
            Matrix&         iy,
      const Tensor7&        scat_i_p,
      const Tensor7&        scat_i_lat,
      const Tensor7&        scat_i_lon,
      const GridPos&        rte_gp_p,
      const GridPos&        rte_gp_lat,
      const GridPos&        rte_gp_lon,
      const Vector&         rte_los,
      const Index&          cloudbox_on,
      const ArrayOfIndex&   cloudbox_limits,
      const Index&          atmosphere_dim,
      const Index&          stokes_dim,
      const Vector&         scat_za_grid,
      const Vector&         scat_aa_grid,
      const Vector&         f_grid )
{
  iy_interp_cloudbox_field( iy, scat_i_p, scat_i_lat, scat_i_lon, rte_gp_p, 
                            rte_gp_lat, rte_gp_lon, rte_los, cloudbox_on, 
                            cloudbox_limits, atmosphere_dim, stokes_dim, 
                            scat_za_grid, scat_aa_grid, f_grid, "cubic" );
}

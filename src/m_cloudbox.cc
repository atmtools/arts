/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
   Claudia Emde     <claudia.emde@dlr.de>
   Cory Davis       <cory.davis@metservice.com>	   
                         
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

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
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


/* Workspace method: Doxygen documentation will be auto-generated */
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
  nlinspace(scat_data_raw[0].f_grid, 1e9, 3.848043e+13 , 5);  
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


/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManually(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
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


/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManuallyAltitude(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
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
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  if( atmosphere_dim == 2 )
    { throw runtime_error( "The cloud box is not defined for 2D." ); }

  // Check keyword arguments
  if( z1 >= z2 )
    throw runtime_error( 
                        "The altitude in *z1* must be smaller than the altitude in *z2*." );
  /* These cases are in fact handled
  if( z1 <= z_field(0, 0, 0) )
    throw runtime_error( "The altitude in *z1* must be larger than the "
                         "first value in *z_field*." );
  if( z2 >= z_field(z_field.npages()-1, 0, 0) )
    throw runtime_error( "The altitude in *z2* must be smaller than the "
                         "last value in *z_field*." );
  */
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


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_i_fieldSetClearsky(Tensor6& doit_i_field,
                             const Tensor7& scat_i_p,
                             const Tensor7& scat_i_lat,
                             const Tensor7& scat_i_lon,
                             const Vector& f_grid,
                             const Index& f_index,
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const ArrayOfIndex& cloudbox_limits,
                             const Index& atmosphere_dim,
                             //Keyword:
                             const Index& all_frequencies
                             )
{
  
  out2 << "  Interpolate boundary clearsky field to obtain the initial field.\n";
  
  // Initial field only needs to be calculated from clearsky field for the 
  // first frequency. For the next frequencies the solution field from the 
  // previous frequencies is used. 
  if(atmosphere_dim == 1)
    {
       if(f_index == 0 || all_frequencies == true){
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
         
         doit_i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1,
                             1,
                             1,
                             N_za,
                             N_aa,
                             N_i);
         
         doit_i_field = 0.;
         
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
                     
                     VectorView target_field = doit_i_field(Range(joker),
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
       else{// no interpolation is required for other frequencies,
         // but the boundary needs to be set correctly.
         doit_i_field(0, 0, 0, Range(joker), Range(joker), Range(joker))=
           scat_i_p(f_index, 0, 0, 0, Range(joker), Range(joker),
                    Range(joker));
           doit_i_field(doit_i_field.nvitrines()-1, 0, 0, Range(joker), 
                        Range(joker), Range(joker))=
             scat_i_p(f_index, 1, 0, 0, Range(joker), Range(joker),
                      Range(joker));
       }
    }
  else if(atmosphere_dim == 3)
    {
      if (all_frequencies == false)
        throw runtime_error("Error in doit_i_fieldSetClearsky: For 3D "
                            "all_frequencies option is not implemented \n");

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
    
 
      //doit_i_field
      doit_i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
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
                        
                          VectorView target_field = doit_i_field(Range(joker),
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
                          
                          VectorView target_field = doit_i_field
                            (p_index, Range(joker), lon_index,
                             za_index, aa_index, i);
                        
                          ConstVectorView source_field = scat_i_lat
                            (f_index, p_index, Range(joker),    
                             lon_index, za_index, aa_index, i);
                        
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
                        
                          VectorView target_field = doit_i_field(p_index,
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


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_i_fieldSetConst(//WS Output:
                        Tensor6& doit_i_field,
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
                        const Vector& doit_i_field_values)
{
  out2 << "  Set initial field to constant values: " << doit_i_field_values << "\n"; 

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
      out3 << "  atm_dim = 1\n"; 
      
    // Define the size of doit_i_field.
    doit_i_field.resize((cloudbox_limits[1] - cloudbox_limits[0])+1, 1, 1,  N_za,
                   1, N_i);
    doit_i_field = 0.;

    // Loop over all zenith angle directions.
    for (Index za_index = 0; za_index < N_za; za_index++)
      {
        for (Index i = 0; i < stokes_dim; i++)
          { 
            //set the value for the upper boundary
            doit_i_field(cloudbox_limits[1]-cloudbox_limits[0], 0, 0, za_index,
                    0, i) = 
          scat_i_p(0, 1, 0, 0, za_index, 0, i);
            //set the value for the lower boundary 
            doit_i_field(0, 0, 0, za_index, 0, i) =  
          scat_i_p(0, 0, 0, 0, za_index, 0, i);
            for (Index scat_p_index = 1; scat_p_index < cloudbox_limits[1] - 
                   cloudbox_limits[0]; scat_p_index++ )
              // The field inside the cloudbox is set to some arbitrary value.
              doit_i_field(scat_p_index, 0, 0, za_index, 0, i) =  doit_i_field_values[i];
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
      

      
      out3 << "atm_dim = 3\n";      
      doit_i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
                     (cloudbox_limits[3]- cloudbox_limits[2])+1,
                     (cloudbox_limits[5]- cloudbox_limits[4])+1,
                     N_za, 
                     N_aa,
                     N_i);
      
      doit_i_field = 0.;
      

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
                            doit_i_field(cloudbox_limits[1]-cloudbox_limits[0], 
                                    lat_index-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_p(0, 1, lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i);
                            //set the value for the lower pressure boundary 
                            doit_i_field(0, lat_index-cloudbox_limits[2],
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
                            doit_i_field(p_index-cloudbox_limits[0], 
                                    cloudbox_limits[3]-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_lat(0, p_index-cloudbox_limits[0],
                                         1, lon_index-cloudbox_limits[4],
                                         za_index, aa_index, i);
                            // second boundary
                            doit_i_field(p_index-cloudbox_limits[0], 0, 
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
                            doit_i_field(p_index-cloudbox_limits[0],
                                    lat_index-cloudbox_limits[2],
                                    cloudbox_limits[5]-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_lon(0, p_index-cloudbox_limits[0],
                                         lat_index-cloudbox_limits[2], 1,
                                         za_index, aa_index, i);
                            // second boundary
                            doit_i_field(p_index-cloudbox_limits[0],  
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
                               doit_i_field(p_index-cloudbox_limits[0],
                                       lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i) =  
                                 doit_i_field_values[i];
                             }
                         }
                     }
                } // stokes loop
            } // aa_grid loop
        } // za_grid loop
       
    } // atmosphere dim = 3
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeInit( //WS Output:
                      ArrayOfSingleScatteringData& scat_data_raw,
                      ArrayOfGField3& pnd_field_raw
                      )
{
  scat_data_raw.reserve(20);
  pnd_field_raw.reserve(20); 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAddAll( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfGField3&  pnd_field_raw,
                 // WS Input(needed for checking the datafiles):
                 const Index& atmosphere_dim,
                 const Vector& f_grid,
                 const Vector& p_grid,
                 const Vector& lat_grid,
                 const Vector& lon_grid,
                 const ArrayOfIndex& cloudbox_limits,
                 // Keywords:
                 const String& filename_scat_data,
                 const String& pnd_field_file)
{
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  
  // Cloudbox limits
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains"
                        "the upper and lower\n"
                        "limit of the cloud for all atmospheric dimensions.\n"
                        "So its length must be 2 x *atmosphere_dim*" ); 
  // Frequency grid
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  
  
  //--- Reading the data ---------------------------------------------------
  ArrayOfString data_files;
  xml_read_from_file(filename_scat_data, data_files);
  scat_data_raw.resize(data_files.nelem());
  
  for (Index i = 0; i<data_files.nelem(); i++)
    {
      
      out2 << "  Read single scattering data\n";
      xml_read_from_file( data_files[i], 
                          scat_data_raw[i]);
      
      chk_single_scattering_data(scat_data_raw[i],
                                 data_files[i], f_grid);  
      
    }
  
  out2 << "  Read particle number density date \n";
  xml_read_from_file(pnd_field_file, pnd_field_raw);
  
  chk_pnd_raw_data(pnd_field_raw,
                   pnd_field_file, atmosphere_dim, p_grid, lat_grid,
                   lon_grid, cloudbox_limits);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAdd( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfGField3&  pnd_field_raw,
                 // WS Input (needed for checking the datafiles):
                 const Index& atmosphere_dim,
                 const Vector& f_grid,
                 const Vector& p_grid,
                 const Vector& lat_grid,
                 const Vector& lon_grid,
                 const ArrayOfIndex& cloudbox_limits,
                 // Keywords:
                 const String& scat_data_file,
                 const String& pnd_field_file
          )
{
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Cloudbox limits
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains"
                        "the upper and lower\n"
                        "limit of the cloud for all atmospheric dimensions.\n"
                        "So its length must be 2 x *atmosphere_dim*" ); 
  // Frequency grid
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  

  //--- Reading the data ---------------------------------------------------

  // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors. 
  SingleScatteringData scat_data;
  scat_data_raw.push_back(scat_data);
  
  GField3 pnd_field_data;
  pnd_field_raw.push_back(pnd_field_data);
  
  out2 << "  Read single scattering data\n";
  xml_read_from_file(scat_data_file, scat_data_raw[scat_data_raw.nelem()-1]);

  chk_single_scattering_data(scat_data_raw[scat_data_raw.nelem()-1],
                             scat_data_file, f_grid);       
  
  out2 << "  Read particle number density field\n";
  if (pnd_field_file.nelem() < 1)
    out1 << "Warning: No pnd_field_file specified. Ignored. \n";
  else
    {
      xml_read_from_file(pnd_field_file,
                         pnd_field_raw[pnd_field_raw.nelem()-1]);
      
      chk_pnd_data(pnd_field_raw[pnd_field_raw.nelem()-1],
                   pnd_field_file, atmosphere_dim, p_grid, lat_grid,
                   lon_grid, cloudbox_limits);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalc(//WS Output:
                   Tensor4& pnd_field,
                   //WS Input
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const ArrayOfGField3& pnd_field_raw,
                   const Index& atmosphere_dim,
                   const ArrayOfIndex& cloudbox_limits
                   )
{
  // Basic checks of input variables
  //
  // Particle number density data
  // 
  if (pnd_field_raw.nelem() == 0)
    throw runtime_error(
                        "No particle number density data given. Please\n"
                        "use WSMs *ParticleTypeInit* and \n"
                        "*ParticleTypeAdd(All)* for reading cloud particle\n"
                        "data.\n"
                        );
  
  // Atmosphere
  if (!((atmosphere_dim == 1) || (atmosphere_dim == 3)))
    throw runtime_error(
                        "*atmosphere_dim* must be 1 or 3, because \n"
                        "scattering is only implemented for 1D and 3D.\n"
                        "DOIT works in 1D and 3D, Monte Carlo only in 3D.\n"
                        );
      
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  const Index Np_cloud = cloudbox_limits[1]-cloudbox_limits[0]+1;
  
  ConstVectorView p_grid_cloud = p_grid[Range(cloudbox_limits[0], Np_cloud)];
      
  //==========================================================================
  if ( atmosphere_dim == 1)
    {
      //Resize variables
      pnd_field.resize(pnd_field_raw.nelem(), Np_cloud, 1, 1);
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
         
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID), p_grid_cloud);
         
          // Interpolation weights:
          Matrix itw(Np_cloud, 2);
          // (2 interpolation weights are required for 1D interpolation)
          interpweights( itw, gp_p);
          // Interpolate:
          interp( pnd_field(i, joker, 0, 0),
                  itw, pnd_field_raw[i](joker, 0, 0), gp_p);
        }
      
    }
  //=========================================================================
  // (CE: Commented atmosphere_dim = 2, becaue scattering is only implemented 
  // in 1D and 3D
  //  else if(atmosphere_dim == 2)
  //     {
  //       //Resize variables
  //       pnd_field.resize(pnd_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
  //                        1);
      
      
  //       // Gridpositions:
  //       ArrayOfGridPos gp_p(p_grid.nelem());
  //       ArrayOfGridPos gp_lat(lat_grid.nelem());
                  
  //       // Interpolate pnd_field. 
  //       // Loop over the particle types:
  //       for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
  //         {
  //           // Calculate grid positions:
  //           p2gridpos(gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID), p_grid);
  //           gridpos(gp_lat, pnd_field_raw[i](0, joker, 0), 
  //                   lat_grid);
      
  //           // Interpolation weights:
  //       Tensor3 itw(p_grid.nelem(), lat_grid.nelem(), 4);
  //       // (8 interpolation weights are required for 3D interpolation)
  //           interpweights( itw, gp_p, gp_lat);
          
  //           // Interpolate:
  //           interp( pnd_field(i, joker, joker, 0),
  //                   itw, pnd_field_raw[i](joker, joker, 0),
  //                   gp_p, gp_lat);
  //         }
  //     }

  //================================================================
  // atmosphere_dim = 3    
  //
  else
    {
      const Index Nlat_cloud = cloudbox_limits[3]-cloudbox_limits[2]+1;
      const Index Nlon_cloud = cloudbox_limits[5]-cloudbox_limits[4]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits[2], Nlat_cloud)];           
      ConstVectorView lon_grid_cloud = 
        lon_grid[Range(cloudbox_limits[4], Nlon_cloud)];
      
      //Resize variables
      pnd_field.resize(pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 
                       Nlon_cloud);
      
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      ArrayOfGridPos gp_lon(Nlon_cloud);
      
      
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                    p_grid_cloud);
          gridpos(gp_lat, pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
                  lat_grid_cloud);
          gridpos(gp_lon, pnd_field_raw[i].get_numeric_grid(GFIELD3_LON_GRID),
                  lon_grid_cloud);
          
          // Interpolation weights:
          Tensor4 itw(Np_cloud, Nlat_cloud, Nlon_cloud, 8);
          // (8 interpolation weights are required for 3D interpolation)
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( pnd_field(i, joker, joker, joker),
                  itw, pnd_field_raw[i], gp_p, gp_lat, gp_lon);
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void DoitCloudboxFieldPut(//WS Output:
               Tensor7&  scat_i_p,
               Tensor7& scat_i_lat,
               Tensor7& scat_i_lon,
               Tensor4& doit_i_field1D_spectrum,
               //WS Input:
               const Tensor6& doit_i_field,
               const Vector& f_grid,
               const Index& f_index,
               const Vector& p_grid,
               const Vector& lat_grid,
               const Vector& lon_grid,
               const Vector& scat_za_grid,
               const Vector& scat_aa_grid,
               const Index& stokes_dim,
               const Index& atmosphere_dim,
               const ArrayOfIndex& cloudbox_limits,
               const Matrix& sensor_pos,
               const Tensor3& z_field
               )
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
  
  // Resize and initialize *doit_i_field_spectra*
  doit_i_field1D_spectrum.resize(N_f, N_p, N_za, stokes_dim); 
  doit_i_field1D_spectrum = 0;

  // Put the doit_i_field at the cloudbox boundary into the interface variable 
  // scat_i_p.
  if(atmosphere_dim == 1)
    {
      bool in_cloudbox = false;
      // Check if sensor inside the cloudbox:
      //loop over all sensor positions
      for (Index i = 0; i < sensor_pos.ncols(); i++)
        {
          // ??? (CE) I think this is worng 
          if(sensor_pos(i, 0) >= z_field(cloudbox_limits[0], 0, 0) &&
             sensor_pos(i, 0) <= z_field(cloudbox_limits[1], 0, 0) )
            {
              in_cloudbox = true;
              out2 << "  Sensor position in cloudbox, store radiation field\n"
                   << "  in cloudbox for all frequencies. \n"; 
            }
        }
      
      // Check size of doit_i_field.
      assert ( is_size( doit_i_field, 
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
              
              //doit_i_field at lower boundary
              scat_i_p(f_index, 0, 0, 0,
                       za, 0, i) = 
                doit_i_field(0, 0, 0, za, 0, i);
              //doit_i_field at upper boundary
              scat_i_p(f_index, 1, 0, 0,
                       za, 0, i) = 
                doit_i_field(cloudbox_limits[1] - cloudbox_limits[0],
                        0, 0, za, 0, i); 

              // If a sensor pos is inside the cloudbox we also need to 
              // define *doit_i_field1D_spectra*
              if( in_cloudbox)
                {
                  doit_i_field1D_spectrum(f_index, joker, za, i) = 
                    doit_i_field(joker, 0, 0, za, 0, i);
                }
              
            }//end stokes_dim
        }//end za loop
      
      
    }//end atmosphere_dim = 1
        
  if(atmosphere_dim == 3)
    {
      // Some sizes relevant for 3D atmosphere
      Index N_lat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      Index N_lon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      
      // Check size of doit_i_field.
      assert ( is_size( doit_i_field, 
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
                  // Put doit_i_field in scat_i_p:
                  //
                  for (Index lat = 0; lat < N_lat; lat++)
                    {
                      for (Index lon = 0; lon < N_lon; lon++)
                        {
                          //doit_i_field at lower boundary
                          scat_i_p(f_index, 0, lat, lon,
                                   za, aa, i) = 
                            doit_i_field(0, lat, lon, za, aa, i);
                          //doit_i_field at upper boundary
                          scat_i_p(f_index, 1, lat, lon,
                                   za, aa, i) = 
                            doit_i_field(cloudbox_limits[1]-cloudbox_limits[0],
                                    lat, lon, za, aa, i);
                        }
                    }
                  // 
                  // Put doit_i_field in scat_i_lat:
                  //
                  for (Index p = 0; p < N_p; p++)
                    {
                      for (Index lon = 0; lon < N_lon; lon++)
                        {
                          //doit_i_field at lower boundary
                          scat_i_lat(f_index, p, 0, lon,
                                     za, aa, i) = 
                            doit_i_field(p, 0, lon, za, aa, i);
                          //doit_i_field at upper boundary
                          scat_i_lat(f_index, p, 1, lon,
                                     za, aa, i) = 
                            doit_i_field(p, cloudbox_limits[3]-
                                    cloudbox_limits[2],
                                    lon, za, aa, i);
                        }
                      //
                      // Put doit_i_field in scat_i_lon:
                      for (Index lat = 0; lat < N_lat; lat++)
                        {
                          //doit_i_field at lower boundary
                          scat_i_lon(f_index, p, lat, 0,
                                     za, aa, i) = 
                            doit_i_field(p, lat, 0, za, aa, i);
                          //doit_i_field at upper boundary
                          scat_i_lon(f_index, p, lat, 1,
                                     za, aa, i) = 
                            doit_i_field(p, lat, cloudbox_limits[5]-
                                    cloudbox_limits[4], za, aa, i);
                        } 
                    }
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void CloudboxGetIncoming(
              Tensor7&        scat_i_p,
              Tensor7&        scat_i_lat,
              Tensor7&        scat_i_lon,
              Index&          cloudbox_on,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         surface_prop_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        const Tensor4&        vmr_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Vector&         scat_za_grid,
        const Vector&         scat_aa_grid )
{
  Index  Nf       = f_grid.nelem();
  Index  Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index  Nza      = scat_za_grid.nelem();
  Index  Ni       = stokes_dim;
  Matrix iy;
  Ppath  ppath;


  //--- Check input ----------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  //if( cloudbox_on == 0  ||  cloudbox_limits.nelem() == 0 )
  // throw runtime_error( "The cloudbox must be activated, and it is not.");
  if( scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180. )
        throw runtime_error(
                     "*scat_za_grid* must include 0 and 180 degrees as endpoints." );
  //--------------------------------------------------------------------------


  // Dummy variable for flag cloudbox_on. It has to be 0 here not to get
  // stuck in an infinite loop (if some propagation path hits the cloud
  // box at some other position.
  cloudbox_on = 0;

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

      //--- Get scat_i_p at lower and upper boundary
      //    (boundary=0: lower, boundary=1: upper)
      for (Index boundary = 0; boundary <= 1; boundary++)
        {
          pos[0] = r_geoid(0,0) + z_field( cloudbox_limits[boundary], 0, 0 );

          for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index ++)
            {
              los[0] =  scat_za_grid[scat_za_index];

              iy_calc_no_jacobian( iy, ppath,
                                   ppath_step_agenda, rte_agenda, 
                                   iy_space_agenda, surface_prop_agenda,
                                   iy_cloudbox_agenda, atmosphere_dim,
                                   p_grid, lat_grid, lon_grid, z_field,
                                   t_field, vmr_field,
                                   r_geoid, z_surface, cloudbox_on,
                                   cloudbox_limits, pos, los, f_grid,
                                   stokes_dim, agenda_verb );

              scat_i_p( joker, boundary, 0, 0, scat_za_index, 0, joker ) = iy;
            }
        }
    }
  

  //--- atmosphere_dim = 3: --------------------------------------------------
  else
    {
      Index Naa = scat_aa_grid.nelem();

      if( scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360. )
        throw runtime_error(
                     "*scat_aa_grid* must include 0 and 360 degrees as endpoints." );

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

      
      //--- Get scat_i_p at lower and upper boundary
      //    (boundary=0: lower, boundary=1: upper)
      for (Index boundary = 0; boundary <= 1; boundary++)
        {
          for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
            {
              for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
                {
                  pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                                   lon_index + cloudbox_limits[4]) 
                    + z_field(cloudbox_limits[boundary],
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
                              iy_calc_no_jacobian( iy, ppath, 
                                                   ppath_step_agenda, 
                                                   rte_agenda, iy_space_agenda,
                                                   surface_prop_agenda, 
                                                   iy_cloudbox_agenda,
                                                   atmosphere_dim, p_grid, 
                                                   lat_grid, lon_grid, z_field,
                                                   t_field, vmr_field,
                                                   r_geoid, z_surface, 
                                                   cloudbox_on,
                                                   cloudbox_limits, 
                                                   pos, los, f_grid,
                                                   stokes_dim, agenda_verb );
                            }

                          scat_i_p( joker, boundary, lat_index, lon_index, 
                                    scat_za_index, scat_aa_index, joker) = iy;
                        }
                    }
                }
            }
        }
      
      //--- Get scat_i_lat (2nd and 3rd boundary)
      for (Index boundary = 0; boundary <= 1; boundary++)
        {
          for (Index p_index = 0; p_index < Np_cloud; p_index++ )
            {
              for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
                {
                  pos[0] = r_geoid(cloudbox_limits[boundary+2],
                                   lon_index + cloudbox_limits[4])
                    + z_field(p_index + cloudbox_limits[0],
                              cloudbox_limits[boundary+2],
                              lon_index + cloudbox_limits[4]);
                  pos[1] = lat_grid[cloudbox_limits[boundary+2]];
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
                              iy_calc_no_jacobian( iy, ppath, 
                                                   ppath_step_agenda, 
                                                   rte_agenda, iy_space_agenda,
                                                   surface_prop_agenda, 
                                                   iy_cloudbox_agenda,
                                                   atmosphere_dim, p_grid, 
                                                   lat_grid, lon_grid, z_field,
                                                   t_field, vmr_field,
                                                   r_geoid, z_surface, 
                                                   cloudbox_on,
                                                   cloudbox_limits, pos, los,
                                                   f_grid, stokes_dim,
                                                   agenda_verb );
                            }

                          scat_i_lat( joker, p_index, boundary, lon_index, 
                                      scat_za_index, scat_aa_index, joker) = iy;
                        }
                    }
                }
            }    
        }

      //--- Get scat_i_lon (1st and 2nd boundary):
      for (Index boundary = 0; boundary <= 1; boundary++)
        {
          for (Index p_index = 0; p_index < Np_cloud; p_index++ )
            {
              for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
                {
                  pos[0] = r_geoid(lat_index + cloudbox_limits[2],
                                   cloudbox_limits[boundary+4]) 
                    + z_field(p_index + cloudbox_limits[0],
                              lat_index + cloudbox_limits[2],
                              cloudbox_limits[boundary+4]);
                  pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
                  pos[2] = lon_grid[cloudbox_limits[boundary+4]];

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
                              iy_calc_no_jacobian( iy, ppath,
                                                   ppath_step_agenda, 
                                                   rte_agenda, iy_space_agenda,
                                                   surface_prop_agenda, 
                                                   iy_cloudbox_agenda,
                                                   atmosphere_dim, p_grid, 
                                                   lat_grid, lon_grid, z_field,
                                                   t_field, vmr_field, r_geoid,
                                                   z_surface, cloudbox_on,
                                                   cloudbox_limits, 
                                                   pos, los, f_grid,
                                                   stokes_dim, agenda_verb );
                            }

                          scat_i_lon( joker, p_index, lat_index, boundary, 
                                      scat_za_index, scat_aa_index, joker) = iy;
                        }
                    }
                }
            }
        } 
    }// End atmosphere_dim = 3.
  cloudbox_on = 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void CloudboxGetIncoming1DAtm(
              Tensor7&        scat_i_p,
              Tensor7&        scat_i_lat,
              Tensor7&        scat_i_lon,
              Index&          cloudbox_on,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         surface_prop_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        const Tensor4&        vmr_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Vector&         scat_za_grid,
        const Vector&         scat_aa_grid )
{
  Index  Nf       = f_grid.nelem();
  Index  Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index  Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Index  Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
  Index  Nza      = scat_za_grid.nelem();
  Index  Naa      = scat_aa_grid.nelem();
  Index  Ni       = stokes_dim;
  Matrix iy;
  Ppath  ppath;


  //--- Check input ----------------------------------------------------------
  if( atmosphere_dim != 3 )
    throw runtime_error( "The atmospheric dimensionality must be 3.");
//   if( cloudbox_on == 0  ||  cloudbox_limits.nelem() == 0 )
//     throw runtime_error( "The cloudbox must be activated, and it is not.");
  if( scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180. )
    throw runtime_error(
                     "*scat_za_grid* must include 0 and 180 degrees as endpoints." );
  if( scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360. )
    throw runtime_error(
                     "*scat_aa_grid* must include 0 and 360 degrees as endpoints." );
  //--------------------------------------------------------------------------

  // Dummy variable for flag cloudbox_on. It has to be 0 here not to get
  // stuck in an infinite loop (if some propagation path hits the cloud
  // box at some other position.
  cloudbox_on = 0;

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

  // These variables are constant for all calculations:
  pos[1] = lat_grid[cloudbox_limits[2]];
  pos[2] = lon_grid[cloudbox_limits[4]];
  los[1] = aa_grid[aa_index];
  
  // Calculate scat_i_p, scat_i_lat, scat_i_lon
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      pos[0] = r_geoid(cloudbox_limits[2], cloudbox_limits[4]) 
        + z_field(cloudbox_limits[0] + p_index, cloudbox_limits[2],
                  cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index++ )
        {
          los[0] = scat_za_grid[scat_za_index];
          
          iy_calc_no_jacobian(iy, ppath, ppath_step_agenda,
                              rte_agenda, iy_space_agenda, surface_prop_agenda, 
                              iy_cloudbox_agenda, atmosphere_dim, p_grid,
                              lat_grid, lon_grid, z_field, t_field, vmr_field,
                              r_geoid, z_surface,
                              cloudbox_on, cloudbox_limits, pos, 
                              los, f_grid, stokes_dim, agenda_verb );
          
          for (Index aa = 0; aa < Naa; aa ++)
            {
              // scat_i_p lower boundary
              if(p_index == 0)
                {
                  for (Index lat = 0; lat < Nlat_cloud; lat ++)
                    {
                      for (Index lon = 0; lon < Nlon_cloud; lon ++)
                        {
                          scat_i_p( joker, 0, lat, lon, scat_za_index, aa,
                                    joker )
                            = iy;
                        }
                    }
                }
              //scat_i_p at upper boundary
              else if (p_index == Np_cloud-1)
                for (Index lat = 0; lat < Nlat_cloud; lat ++)
                  {
                    for (Index lon = 0; lon < Nlon_cloud; lon ++)
                      {
                        scat_i_p( joker, 1, lat, lon, scat_za_index, aa,
                                  joker )
                          = iy;
                      }
                  }
              
              // scat_i_lat (both boundaries)
              for (Index lat = 0; lat < 2; lat ++) 
                {
                  for (Index lon = 0; lon < Nlon_cloud; lon ++)
                    {
                      scat_i_lat( joker, p_index, lat, lon, 
                                  scat_za_index, aa, joker)
                        = iy;
                    }
                }
              
              // scat_i_lon (both boundaries)
              for (Index lat = 0; lat < Nlat_cloud; lat ++) 
                {
                  for (Index lon = 0; lon < 2; lon ++)
                    {
                      scat_i_lon( joker, p_index, lat, lon, 
                                  scat_za_index, aa, joker)
                        = iy;
                    }
                }
            }
        }
    }
  cloudbox_on = 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void iyInterpCloudboxField(
            Matrix&         iy,
      const Tensor7&        scat_i_p,
      const Tensor7&        scat_i_lat,
      const Tensor7&        scat_i_lon,
      const Tensor4&        doit_i_field1D_spectrum,
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
      const Vector&         f_grid)
{
  iy_interp_cloudbox_field( iy, scat_i_p, scat_i_lat, scat_i_lon, 
                            doit_i_field1D_spectrum, rte_gp_p, 
                            rte_gp_lat, rte_gp_lon, rte_los, cloudbox_on, 
                            cloudbox_limits, atmosphere_dim, stokes_dim, 
                            scat_za_grid, scat_aa_grid, f_grid, "linear" );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void iyInterpPolyCloudboxField(
            Matrix&         iy,
      const Tensor7&        scat_i_p,
      const Tensor7&        scat_i_lat,
      const Tensor7&        scat_i_lon,
      const Tensor4&        doit_i_field1D_spectrum,      
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
  iy_interp_cloudbox_field( iy, scat_i_p, scat_i_lat, scat_i_lon, 
                            doit_i_field1D_spectrum, rte_gp_p, 
                            rte_gp_lat, rte_gp_lon, rte_los, cloudbox_on, 
                            cloudbox_limits, atmosphere_dim, stokes_dim, 
                            scat_za_grid, scat_aa_grid, f_grid, "polynomial" );
}


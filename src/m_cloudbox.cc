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
  iy_cloudbox_agenda = Agenda();
  iy_cloudbox_agenda.set_name( "iy_cloudbox_agenda" );
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
void ParticleTypeInit( //WS Output:
                      ArrayOfSingleScatteringData& scat_data_raw,
                      ArrayOfGriddedField3& pnd_field_raw
                      )
{
  scat_data_raw.reserve(20);
  pnd_field_raw.reserve(20); 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAddAll( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfGriddedField3&  pnd_field_raw,
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
                 ArrayOfGriddedField3&  pnd_field_raw,
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
  
  GriddedField3 pnd_field_data;
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
                   const ArrayOfGriddedField3& pnd_field_raw,
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
      pnd_field.resize(pnd_field_raw.nelem(), Np_cloud, 1, 1 );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
         
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos( gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID), 
                     p_grid_cloud );

          // Interpolation weights:
          Matrix itw(Np_cloud, 2);

          interpweights( itw, gp_p);
          // Interpolate:
          interp( pnd_field(i,joker,0,0), itw, 
                  pnd_field_raw[i](joker,0,0), gp_p );
        }
    }

  else if(atmosphere_dim == 2)
    {
      const Index Nlat_cloud = cloudbox_limits[3]-cloudbox_limits[2]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits[2],Nlat_cloud)];           
      
      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 1 );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos( gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                     p_grid_cloud);
          gridpos( gp_lat, pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
                   lat_grid_cloud);
          
          // Interpolation weights:
          Tensor3 itw( Np_cloud, Nlat_cloud, 4 );
          interpweights( itw, gp_p, gp_lat );
          
          // Interpolate:
          interp( pnd_field(i,joker,joker,0), itw, 
                  pnd_field_raw[i](joker,joker,0), gp_p, gp_lat );
        }
    }
  else
    {
      const Index Nlat_cloud = cloudbox_limits[3]-cloudbox_limits[2]+1;
      const Index Nlon_cloud = cloudbox_limits[5]-cloudbox_limits[4]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits[2],Nlat_cloud)];           
      ConstVectorView lon_grid_cloud = 
        lon_grid[Range(cloudbox_limits[4],Nlon_cloud)];
      
      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 
                        Nlon_cloud );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      ArrayOfGridPos gp_lon(Nlon_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos( gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                     p_grid_cloud);
          gridpos( gp_lat, pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
                   lat_grid_cloud);
          gridpos( gp_lon, pnd_field_raw[i].get_numeric_grid(GFIELD3_LON_GRID),
                   lon_grid_cloud);
          
          // Interpolation weights:
          Tensor4 itw( Np_cloud, Nlat_cloud, Nlon_cloud, 8 );
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( pnd_field(i,joker,joker,joker), itw, 
                  pnd_field_raw[i], gp_p, gp_lat, gp_lon );
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldZero(//WS Output:
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


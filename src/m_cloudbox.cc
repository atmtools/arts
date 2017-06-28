/* Copyright (C) 2002-2012
   Patrick Eriksson <patrick.eriksson@chalmers.se>
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
#include <cmath>
#include <stdexcept>
#include <cstdlib>

#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "file.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric LAT_LON_MIN;
extern const String SCATSPECIES_MAINTAG;
extern const Numeric DENSITY_OF_ICE;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxOff (
         Index&                       cloudbox_on,
         ArrayOfIndex&                cloudbox_limits,
         Agenda&                      iy_cloudbox_agenda,
         Tensor4&                     pnd_field,
         ArrayOfArrayOfSingleScatteringData& scat_data,
         Matrix&                      particle_masses,
   const Verbosity&)
{
  cloudbox_on = 0;
  cloudbox_limits.resize ( 0 );
  iy_cloudbox_agenda = Agenda();
  iy_cloudbox_agenda.set_name ( "iy_cloudbox_agenda" );
  pnd_field.resize(0,0,0,0);
  scat_data.resize(0);
  particle_masses.resize(0,0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetAutomatically (// WS Output:
                               //Workspace& /* ws */,
                               Index&          cloudbox_on,
                               ArrayOfIndex&   cloudbox_limits,
                               //Agenda&  iy_cloudbox_agenda,
                               // WS Input:
                               const Index&    atmosphere_dim,
                               const Vector&   p_grid,
                               const Vector&   lat_grid,
                               const Vector&   lon_grid,
                               const Tensor4&  scat_species_mass_density_field,
                               const Tensor4&  scat_species_mass_flux_field,
                               const Tensor4&  scat_species_number_density_field,
                               const Tensor4&  scat_species_mean_mass_field,
                               // Control Parameters
                               const Numeric&  cloudbox_margin,
                               const Verbosity& verbosity)
{
  // Check existing WSV
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  // includes p_grid chk_if_decresasing
  chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid ); 
  // Set cloudbox_on
  cloudbox_on = 1;

  if ( atmosphere_dim > 1 )
    {
      ostringstream os;
      os << "cloudboxSetAutomatically not yet available for 2D and 3D cases.";
      throw runtime_error( os.str() );
    }

  // Allocate cloudbox_limits
  cloudbox_limits.resize ( atmosphere_dim*2 );

  // Variables
  Index p1;
  if ( cloudbox_margin == -1 )
    {
      cloudbox_limits[0] = 0;
      p1 = 0;
    }
  else p1 = scat_species_mass_density_field.npages()-1;
  Index p2 = 0;

  // OLE: Commented out until code that uses it at the end of this function is commented back in
//  if ( atmosphere_dim > 1 )
//    {
//      Index lat1 = scat_species_mass_density_field.nrows()-1;
//      Index lat2 = 0;
//    }
//  if ( atmosphere_dim > 2 )
//    {
//      Index lon1 = scat_species_mass_density_field.ncols()-1;
//      Index lon2 = 0;
//    }

  bool not_empty_any=false;
  bool not_empty_md=true;
  bool not_empty_mf=true;
  bool not_empty_nd=true;
  bool not_empty_mm=true;

  Index nss=0;

  if (scat_species_mass_density_field.empty())
    not_empty_md = 0;
  else
    nss = scat_species_mass_density_field.nbooks();

  if (scat_species_mass_flux_field.empty())
    not_empty_mf = 0;
  else if (nss!=0)
    {
      if (nss!=scat_species_mass_flux_field.nbooks())
        {
          ostringstream os;
          os << "Inconsistent number of scattering elements in\n"
             << "scat_species_mass_flux_field compared to other\n"
             << "scat.species fields.";
          throw runtime_error( os.str() );
        }
     }
  else
    nss = scat_species_mass_flux_field.nbooks();

  if (scat_species_number_density_field.empty())
    not_empty_nd = 0;
  else if (nss!=0)
    {
      if (nss!=scat_species_number_density_field.nbooks())
        {
          ostringstream os;
          os << "Inconsistent number of scattering elements in\n"
             << "scat_species_number_density_field compared to other\n"
             << "scat.species fields.";
          throw runtime_error( os.str() );
        }
     }
  else
    nss = scat_species_mass_flux_field.nbooks();

  if (scat_species_mean_mass_field.empty())
    not_empty_mm = 0;
  else if (nss!=0)
    {
      if (nss!=scat_species_mean_mass_field.nbooks())
        {
          ostringstream os;
          os << "Inconsistent number of scattering elements in\n"
             << "scat_species_mean_mass_field compared to other\n"
             << "scat.species fields.";
          throw runtime_error( os.str() );
        }
     }
  else
    nss = scat_species_mean_mass_field.nbooks();


  //--------- Start loop over scattering species ------------------------------
  for ( Index l=0; l<nss; l++ )
  {
    //cout << "for scatt species #" << l << ":\n";
    bool not_empty;

    //not_empty is set to true, if a single value of scat_species_XX_field
    //is unequal zero (and not NaN), i.e. if we actually have some amount of
    //these scattering species in the atmosphere.
    if (not_empty_md)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_mass_density_field ( l, joker, joker, joker ),
                               "scat_species_mass_density_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      //if particles found, enter detailed search
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_mass_density_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in mass density field: " << not_empty << "\n";
    }

    if (not_empty_mf)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_mass_flux_field ( l, joker, joker, joker ),
                               "scat_species_mass_flux_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_mass_flux_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in mass flux field: " << not_empty << "\n";
    }

    if (not_empty_nd)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_number_density_field ( l, joker, joker, joker ),
                               "scat_species_number_density_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_number_density_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in number density field: " << not_empty << "\n";
    }

    if (not_empty_mm)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_mean_mass_field ( l, joker, joker, joker ),
                               "scat_species_mean_mass_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_mean_mass_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in mean mass field: " << not_empty << "\n";
    }
    //cout << "particles in any field: " << not_empty_any << "\n";
  }

  // decrease lower cb limit by one to ensure that linear interpolation of 
  // particle number densities is possible.
  Index p0 = 0; //only for the use of function *max*
  p1 = max(p1-1, p0);

  Numeric p_margin1;

  // alter lower cloudbox_limit by cloudbox_margin, using barometric
  // height formula
  p_margin1 = barometric_heightformula ( p_grid[p1], cloudbox_margin );
  Index k = 0;
  while ( p_grid[k+1] >= p_margin1 && k+1 < p_grid.nelem() ) k++;
  cloudbox_limits[0]= k;

  // increase upper cb limit by one to ensure that linear interpolation of 
  // particle number densities is possible.
  p2 = min(p2+1, scat_species_mass_density_field.npages()-1);
  // set upper cloudbox_limit
  // if cloudbox reaches to the upper most pressure level
  if ( p2 >= scat_species_mass_density_field.npages()-1)
  {
    CREATE_OUT2;
    out2<<"The cloud reaches to TOA!\n"
    <<"Check scat_species_mass_density_field data, if realistic!\n";
  }
  cloudbox_limits[1] = p2;

  //out0<<"\n"<<p2<<"\n"<<p_grid[p2]<<"\n";

  // check if all selected scattering species fields are zero at each level,
  // than switch cloudbox off, skipping scattering calculations
  if ( !not_empty_any )
  {
    CREATE_OUT0;
    //cloudboxOff ( cloudbox_on, cloudbox_limits, iy_cloudbox_agenda );
    cloudbox_on = 0;
    out0<<"Cloudbox is switched off!\n";

    return;
  }


  // assert keyword arguments

  // The pressure in *p1* must be bigger than the pressure in *p2*.
  assert ( p_grid[p1] > p_grid[p2] );
  // The pressure in *p1* must be larger than the last value in *p_grid*.
  assert ( p_grid[p1] > p_grid[p_grid.nelem()-1] );
  // The pressure in *p2* must be smaller than the first value in *p_grid*."
  assert ( p_grid[p2] < p_grid[0] );

  /*
  if ( atmosphere_dim >= 2 )
  {
    // The latitude in *lat2* must be bigger than the latitude in *lat1*.
    assert ( lat_grid[lat2] > lat_grid[lat1] );
    // The latitude in *lat1* must be >= the second value in *lat_grid*.
    assert ( lat_grid[lat1] >= lat_grid[1] );
    // The latitude in *lat2* must be <= the next to last value in *lat_grid*.
    assert ( lat_grid[lat2] <= lat_grid[lat_grid.nelem()-2] );
  }
  if ( atmosphere_dim == 3 )
  {
    // The longitude in *lon2* must be bigger than the longitude in *lon1*.
    assert ( lon_grid[lon2] > lon_grid[lon1] );
    // The longitude in *lon1* must be >= the second value in *lon_grid*.
    assert ( lon_grid[lon1] >= lon_grid[1] );
    // The longitude in *lon2* must be <= the next to last value in *lon_grid*.
    assert ( lon_grid[lon2] <= lon_grid[lon_grid.nelem()-2] );
  }
  */
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetFullAtm(//WS Output
                       Index& cloudbox_on,
                       ArrayOfIndex& cloudbox_limits,
                       // WS Input
                       const Index& atmosphere_dim,
                       const Vector& p_grid,
                       const Vector& lat_grid,
                       const Vector& lon_grid,
                       const Verbosity&)
{
  cloudbox_on = 1;
  cloudbox_limits.resize(2*atmosphere_dim); 

  cloudbox_limits[0] = 0;
  cloudbox_limits[1] = p_grid.nelem()-1;

  if( atmosphere_dim > 1 )
    {
      Index last_lat = lat_grid.nelem()-1;

      // find minimum lat_grid point i with lat_grid[i]-lat_grid[0]>=LAT_LON_MIN
      Index i=1;
      while( (i<last_lat-1) && (lat_grid[i]-lat_grid[0] < LAT_LON_MIN) )
        i++;
      if( i==last_lat-1 )
        {
          ostringstream os;
          os << "Can not define lower latitude edge of cloudbox:\n"
             << "Extend of atmosphere too small. Distance to minimum latitude\n"
             << "has to be at least " << LAT_LON_MIN << "deg, but only "
             << lat_grid[i-1]-lat_grid[0] << " available here.";
          throw runtime_error( os.str() );
        }
      cloudbox_limits[2] = i;

      // find maximum lat_grid point j with lat_grid[-1]-lat_grid[j]>=LAT_LON_MIN
      // and j>i
      Index j=last_lat-1;
      while( (j>i) && (lat_grid[last_lat]-lat_grid[j] < LAT_LON_MIN) )
        j--;
      if( j==i )
        {
          ostringstream os;
          os << "Can not define upper latitude edge of cloudbox:\n"
             << "Extend of atmosphere too small. Distance to maximum latitude\n"
             << "has to be at least " << LAT_LON_MIN << "deg, but only "
             << lat_grid[last_lat]-lat_grid[j+1] << " available here.";
          throw runtime_error( os.str() );
        }
      cloudbox_limits[3] = j;
    }

  if( atmosphere_dim > 2 )
    {
      const Numeric latmax = max( abs(lat_grid[cloudbox_limits[2]]),
                                  abs(lat_grid[cloudbox_limits[3]]) );
      const Numeric lfac = 1 / cos( DEG2RAD*latmax );
      Index last_lon = lon_grid.nelem()-1;

      // find minimum lon_grid point i with lon_grid[i]-lon_grid[0]>=LAT_LON_MIN/lfac
      Index i=1;
      while( (i<last_lon-1) && (lon_grid[i]-lon_grid[0] < LAT_LON_MIN/lfac) )
        i++;
      if( i==last_lon-1 )
        {
          ostringstream os;
          os << "Can not define lower longitude edge of cloudbox:\n"
             << "Extend of atmosphere too small. Distance to minimum longitude\n"
             << "has to be at least " << LAT_LON_MIN/lfac << "deg, but only "
             << lon_grid[i-1]-lon_grid[0] << " available here.";
          throw runtime_error( os.str() );
        }
      cloudbox_limits[4] = i;

      // find maximum lon_grid point j with lon_grid[-1]-lon_grid[j]>=LAT_LON_MIN/lfac
      // and j>i
      Index j=last_lon-1;
      while( (j>i) && (lon_grid[last_lon]-lon_grid[j] < LAT_LON_MIN/lfac) )
        j--;
      if( j==i )
        {
          ostringstream os;
          os << "Can not define upper longitude edge of cloudbox:\n"
             << "Extend of atmosphere too small. Distance to maximum longitude\n"
             << "has to be at least " << LAT_LON_MIN/lfac << "deg, but only "
             << lon_grid[last_lon]-lon_grid[j+1] << " available here.";
          throw runtime_error( os.str() );
        }
      cloudbox_limits[5] = j;
    }
}
  

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManually(// WS Output:
                         Index&         cloudbox_on,
                         ArrayOfIndex&  cloudbox_limits,
                         // WS Input:
                         const Index&   atmosphere_dim,
                         const Vector&  p_grid,
                         const Vector&  lat_grid,
                         const Vector&  lon_grid,
                         // Control Parameters:
                         const Numeric& p1,
                         const Numeric& p2,
                         const Numeric& lat1,
                         const Numeric& lat2,
                         const Numeric& lon1,
                         const Numeric& lon2,
                         const Verbosity&)
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Check keyword arguments
  if( p1 <= p2 )
    throw runtime_error( "The pressure in *p1* must be bigger than the "
                         "pressure in *p2*." );
  if( p1 <= p_grid[p_grid.nelem()-1] )
    throw runtime_error( "The pressure in *p1* must be larger than the "
                         "last value in *p_grid*." );
  if( p2 >= p_grid[0] )
    throw runtime_error( "The pressure in *p2* must be smaller than the "
                         "first value in *p_grid*." );
  if( atmosphere_dim >= 2 )
    {
      if( lat2 <= lat1 )
        throw runtime_error( "The latitude in *lat2* must be bigger than the "
                             "latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
        throw runtime_error( "The latitude in *lat1* must be >= the "
                             "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
        throw runtime_error( "The latitude in *lat2* must be <= the "
                             "next to last value in *lat_grid*." );
    }
  if( atmosphere_dim == 3 )
    {
      if( lon2 <= lon1 )
        throw runtime_error( "The longitude in *lon2* must be bigger than the "
                             "longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
        throw runtime_error( "The longitude in *lon1* must be >= the "
                             "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
        throw runtime_error( "The longitude in *lon2* must be <= the "
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
void cloudboxSetManuallyAltitude(// WS Output:
                                 Index&         cloudbox_on,
                                 ArrayOfIndex&  cloudbox_limits,
                                 // WS Input:
                                 const Index&   atmosphere_dim,
                                 const Tensor3& z_field,
                                 const Vector&  lat_grid,
                                 const Vector&  lon_grid,
                                 // Control Parameters:
                                 const Numeric& z1,
                                 const Numeric& z2,
                                 const Numeric& lat1,
                                 const Numeric& lat2,
                                 const Numeric& lon1,
                                 const Numeric& lon2,
                                 const Verbosity&)
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Check keyword arguments
  if( z1 >= z2 )
    throw runtime_error( "The altitude in *z1* must be smaller than the "
                         "altitude in *z2*." );
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
        throw runtime_error( "The latitude in *lat2* must be bigger than the "
                             " latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
        throw runtime_error( "The latitude in *lat1* must be >= the "
                             "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
        throw runtime_error( "The latitude in *lat2* must be <= the "
                             "next to last value in *lat_grid*." );
      if( lon2 <= lon1 )
        throw runtime_error( "The longitude in *lon2* must be bigger than the "
                             "longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
        throw runtime_error( "The longitude in *lon1* must be >= the "
                             "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
        throw runtime_error( "The longitude in *lon2* must be <= the "
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
void scat_species_fieldCleanup (//WS Output:
                          Tensor4& scat_species_field_out,
                          //WS Input:
                          const Tensor4& scat_species_field_in,
                          const Numeric& threshold,
                          const Verbosity&)
{
  // First make a copy of scat_species_field_in as it might be identical with
  // scat_species_field_out. And resize scat_species_field_out appropriately.
  Tensor4 infield = scat_species_field_in;

  // Check that scat_species_field contains realistic values
  //(values smaller than the threshold will be set to 0)
  for ( Index i=0; i<scat_species_field_in.nbooks(); i++ )
  {
    for ( Index j=0; j<scat_species_field_in.npages(); j++ )
    {
      for ( Index k=0; k<scat_species_field_in.nrows(); k++ )
      {
        for ( Index l=0; l<scat_species_field_in.ncols(); l++ )
        {
          if ( scat_species_field_in ( i,j,k,l ) < threshold ) 
          {
            infield ( i,j,k,l ) = 0.0;
          }
        }
      }
    }
  }
  scat_species_field_out = infield;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_speciesSet (// WS Generic Output:
                      ArrayOfString&  scat_species,
                      // Control Parameters:
                      const ArrayOfString& scatspecies_tags,
                      const String& delim,
                      const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  scat_species.resize ( scatspecies_tags.nelem() );
  //assign input strings to scat_species
  scat_species = scatspecies_tags;

  chk_scat_species (scat_species, delim);

  // Print list of scattering species settings to the most verbose output stream:
  out3 << "  Defined scattering species settings: ";
  for ( Index i=0; i<scat_species.nelem(); ++i )
  {
    out3 << "\n  " << i << ": "<<scat_species[i];

  }
  out3 << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesInit (//WS Output:
                       ArrayOfString& scat_species,
                       ArrayOfArrayOfSingleScatteringData& scat_data,
                       ArrayOfArrayOfScatteringMetaData& scat_meta,
                       ArrayOfGriddedField3& pnd_field_raw,
                       const Verbosity&)
{
  scat_species.resize ( 0 );
  scat_data.resize(0);
  scat_meta.resize(0);
  pnd_field_raw.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsPndAndScatAdd( //WS Output:
                     ArrayOfArrayOfSingleScatteringData& scat_data,
                     ArrayOfGriddedField3&  pnd_field_raw,
                     // WS Input (needed for checking the datafiles):
                     const Index& atmosphere_dim,
                     // Keywords:
                     const ArrayOfString& scat_data_files,
                     const ArrayOfString& pnd_field_files,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  //chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );


  //--- Reading the data ---------------------------------------------------

  if ( scat_data_files.nelem()!=pnd_field_files.nelem() )
    {
      ostringstream os;
      os << "Number of elements in scat_data and pnd_field filelists is"
         << "inconsistent.";
      throw runtime_error ( os.str() );
    }

  Index last_species = scat_data.nelem()-1;
  if (last_species == -1)
  {
      scat_data.resize(1);
      last_species = 0;
  }

  // Create empty dummy elements to append to *scat_data* and *pnd_field_raw*.
  SingleScatteringData scat_data_single;
  GriddedField3 pnd_field_data;

  for ( Index i = 0; i<scat_data_files.nelem(); i++ )
    {
      // Append *scat_data* and *pnd_field_raw* with empty Arrays of Tensors.
      scat_data[last_species].push_back(scat_data_single);
      pnd_field_raw.push_back(pnd_field_data);
  
      out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
      xml_read_from_file(scat_data_files[i],
                         scat_data[last_species][scat_data[last_species].nelem()-1],
                         verbosity);
  
      out2 << "  Read particle number density field\n";
      if (pnd_field_files[i].nelem() < 1)
      {
        CREATE_OUT1;
        out1 << "Warning: No pnd_field_file specified. Ignored here,\n"
             << "but user HAS TO add that later on!. \n";
      }
      else
      {
        xml_read_from_file(pnd_field_files[i],
                           pnd_field_raw[pnd_field_raw.nelem()-1], verbosity);
      
        chk_pnd_data(pnd_field_raw[pnd_field_raw.nelem()-1],
                     pnd_field_files[i], atmosphere_dim, verbosity);
      }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesPndAndScatAdd (//WS Output:
                         ArrayOfArrayOfSingleScatteringData& scat_data,
                         ArrayOfGriddedField3&  pnd_field_raw,
                         // WS Input(needed for checking the datafiles):
                         const Index& atmosphere_dim,
                         // Keywords:
                         const ArrayOfString& scat_data_files,
                         const String& pnd_fieldarray_file,
                         const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  //chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid );


  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  arr_ssd.resize ( scat_data_files.nelem() );

  for ( Index i = 0; i<scat_data_files.nelem(); i++ )
  {

    out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
    xml_read_from_file ( scat_data_files[i], arr_ssd[i], verbosity );

  }

  // append as new scattering species
  if ( scat_data.nelem()==0 )
    {
      scat_data.resize(1);
      scat_data[0] = arr_ssd;
    }
  else
    scat_data.push_back(arr_ssd);

  out2 << "  Read particle number density data \n";
  ArrayOfGriddedField3  pnd_tmp;
  xml_read_from_file ( pnd_fieldarray_file, pnd_tmp, verbosity );

  chk_pnd_raw_data ( pnd_tmp,
                     pnd_fieldarray_file, atmosphere_dim, verbosity);

  // append to pnd_field_raw
  for ( Index i = 0; i<pnd_tmp.nelem(); ++i)
    pnd_field_raw.push_back(pnd_tmp[i]);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsToabs_speciesAdd( //WS Output:
                     ArrayOfArrayOfSingleScatteringData& scat_data,
                     ArrayOfGriddedField3& vmr_field_raw,
                     ArrayOfArrayOfSpeciesTag& abs_species,
                     Index& propmat_clearsky_agenda_checked,
                     Index& abs_xsec_agenda_checked,
                     // WS Input (needed for checking the datafiles):
                     const Index& atmosphere_dim,
                     const Vector& f_grid,
                     // Keywords:
                     const ArrayOfString& scat_data_files,
                     const ArrayOfString& pnd_field_files,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  //chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Frequency grid
  if( f_grid.empty() )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  

  //--- Reading the data ---------------------------------------------------

  if ( scat_data_files.nelem()!=pnd_field_files.nelem() )
    {
      ostringstream os;
      os << "Number of elements in scat_data and pnd_field filelists is"
         << "inconsistent.";
      throw runtime_error ( os.str() );
    }

  Index last_species = scat_data.nelem()-1;
  if (last_species == -1)
  {
      scat_data.resize(1);
      last_species = 0;
  }

  // Create empty dummy elements to append to *scat_data* and *pnd_field_raw*.
  SingleScatteringData scat_data_single;
  GriddedField3 pnd_field_data;
  ArrayOfString species(1);
  species[0] = "particles";

  for ( Index i = 0; i<scat_data_files.nelem(); i++ )
    {
      // Append *scat_data* and *pnd_field_raw* with empty Arrays of Tensors.
      scat_data[last_species].push_back(scat_data_single);
      vmr_field_raw.push_back(pnd_field_data);

      out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
      xml_read_from_file(scat_data_files[i],
                         scat_data[last_species][scat_data[last_species].nelem()-1],
                         verbosity);

      out2 << "  Check single scattering properties\n";
      chk_scat_data_fgrid(scat_data[last_species][scat_data[last_species].nelem()-1],
                          f_grid, "scat_data_single.f_grid to f_grid");

  
      out2 << "  Read particle number density field\n";
      if (pnd_field_files[i].nelem() < 1)
      {
        CREATE_OUT1;
        out1 << "Warning: No pnd_field_file specified. Ignored here,\n"
             << "but user HAS TO add that later on!\n";
      }
      else
      {
        try
        {
          xml_read_from_file(pnd_field_files[i], vmr_field_raw[vmr_field_raw.nelem()-1], verbosity);
        }
        catch(...)
        {
          ArrayOfGriddedField3 tmp;
          try
          {
            xml_read_from_file(pnd_field_files[i], tmp, verbosity);
            if(tmp.nelem() == 1)
            {
              vmr_field_raw[vmr_field_raw.nelem()-1] = tmp[0];
            }
            else 
            {
              std::ostringstream os;
              os << "The file " << pnd_field_files[i] << "\n" <<
              "is neither GriddedField3 nor a 1-long ArrayOfGriddedField3.\n";
              throw std::runtime_error(os.str()); 
            }
          }
          catch(...)
          {
            std::ostringstream os;
            os << "The file " << pnd_field_files[i] << " does not exist or\n" <<
            "its type is neither GriddedField3 nor a 1-long ArrayOfGriddedField3.\n";
            throw std::runtime_error(os.str());
            
          }
          
        }
      
        chk_pnd_data(vmr_field_raw[vmr_field_raw.nelem()-1],
                     pnd_field_files[i], atmosphere_dim, verbosity);
      }

      out2 << "  Append 'particle' field to abs_species\n";
      abs_speciesAdd( abs_species,
                      propmat_clearsky_agenda_checked, abs_xsec_agenda_checked,
                      species, verbosity );
    }
  scat_dataCheck( scat_data, "sane", 1e-2, verbosity );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesScatAndMetaRead (//WS Output:
                                 ArrayOfArrayOfSingleScatteringData& scat_data,
                                 ArrayOfArrayOfScatteringMetaData& scat_meta,
                                 // Keywords:
                                 const ArrayOfString& scat_data_files,
                                 const Verbosity& verbosity)
{
  CREATE_OUT2;

  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  ArrayOfScatteringMetaData arr_smd;

  arr_ssd.resize ( scat_data_files.nelem() );
  arr_smd.resize ( scat_data_files.nelem() );

  Index meta_naming_conv=0;

  for ( Index i = 0; i<scat_data_files.nelem(); i++ )
    {
      out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
      xml_read_from_file ( scat_data_files[i], arr_ssd[i], verbosity );

      // make meta data name from scat data name
      ArrayOfString strarr;
      String scat_meta_file;

      if( i==0 )
        {
          try
            {
              scat_data_files[i].split ( strarr, ".xml" );
              scat_meta_file = strarr[0]+".meta.xml";

              out2 << "  Read scattering meta data\n";
              xml_read_from_file ( scat_meta_file, arr_smd[i], verbosity );

              meta_naming_conv = 1;
            }
          catch (runtime_error e1)
            {
              try
                {
                  scat_data_files[i].split ( strarr, "scat_data" );
                  scat_meta_file = strarr[0]+"scat_meta"+strarr[1];

                  out2 << "  Read scattering meta data\n";
                  xml_read_from_file ( scat_meta_file, arr_smd[i], verbosity );

                  meta_naming_conv = 2;
                 }
              catch (runtime_error e2)
                {
                  ostringstream os;
                  os << "No meta data file following one of the allowed naming "
                     << "conventions was found.\n"
                     << "Allowed are "
                     << "*.meta.xml from *.xml and "
                     << "*scat_meta* from *scat_data*";
                  throw runtime_error(os.str());
                }
            }
        }
      else
        {
          if( meta_naming_conv==1 )
            {
              scat_data_files[i].split ( strarr, ".xml" );
              scat_meta_file = strarr[0]+".meta.xml";

              out2 << "  Read scattering meta data\n";
              xml_read_from_file ( scat_meta_file, arr_smd[i], verbosity );
            }
          else
            {
              scat_data_files[i].split ( strarr, "scat_data" );
              scat_meta_file = strarr[0]+"scat_meta"+strarr[1];

              out2 << "  Read scattering meta data\n";
              xml_read_from_file ( scat_meta_file, arr_smd[i], verbosity );
            }
        }

      //FIXME: currently nothing is done in chk_scattering_meta_data!
      chk_scattering_meta_data ( arr_smd[i],
                                 scat_meta_file, verbosity );
            
    }

  // check if arrays have same size
  chk_scattering_data ( arr_ssd,
                        arr_smd, verbosity );

  // append as new scattering species
  scat_data.push_back(arr_ssd);
  scat_meta.push_back(arr_smd);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsSelect (//WS Output:
                         ArrayOfArrayOfSingleScatteringData& scat_data,
                         ArrayOfArrayOfScatteringMetaData& scat_meta,
                         // WS Input:
                         const ArrayOfString& scat_species,
                         const String& species,
                         const String& sizeparam,
                         const Numeric& sizemin,
                         const Numeric& sizemax,
                         const Numeric& tolerance,
                         const String& delim,
                         const Verbosity& )
{ 
  // first check that sizes of scat_species and scat_data/scat_meta agree
  Index nspecies = scat_species.nelem();
  if ( nspecies != scat_data.nelem() || nspecies != scat_meta.nelem() )
    {
      ostringstream os;
      os << "Number of scattering species specified by scat_species does\n"
         << "not agree with number of scattering species in\n"
         << "scat_data or scat_meta:\n"
         << "scat_species has " << nspecies << " entries, while scat_data has "
         << scat_data.nelem() << " and scat_meta has " << scat_meta.nelem()
         << ".";
      throw runtime_error ( os.str() );
    }

  // create temporary containers for selected elements
  ArrayOfSingleScatteringData scat_data_tmp;
  ArrayOfScatteringMetaData scat_meta_tmp;

  String partfield_name;
  //find the species to handle: compare 'species' to 'partfield' part of
  //scat_species tags
  Index i_ss=-1;
  for ( Index i=0; i<scat_species.nelem(); i++ )
    {
      parse_partfield_name( partfield_name, scat_species[i], delim);
      if ( partfield_name==species )
        i_ss=i;
    }
  if ( i_ss<0 )
    {
      ostringstream os;
      os << "Scattering species " << species << " not found among scat_species.";
      throw runtime_error ( os.str() );
    }

  // choosing the specified SingleScatteringData and ScatteringMetaData
  if ( sizeparam=="diameter_max" )
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
    {
      // scattering element diameter is extracted from the
      // scattering element's meta data and checked whether it's within size
      // selected range (sizemax < 0 check follows from wildcard usage and
      // means consider all sizes on the upper end)
      if ( scat_meta[i_ss][i_se].diameter_max > sizemin-sizemin*tolerance &&
           ( sizemax+sizemax*tolerance > scat_meta[i_ss][i_se].diameter_max ||
             sizemax < 0. ) )
        {
          // copy selected scattering element to temp arrays
          scat_data_tmp.push_back(scat_data[i_ss][i_se]);
          scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
        }
    }
  else if ( sizeparam=="diameter_volume_equ" )
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
    {
      if ( scat_meta[i_ss][i_se].diameter_volume_equ
             > sizemin-sizemin*tolerance &&
           ( sizemax+sizemax*tolerance
               > scat_meta[i_ss][i_se].diameter_volume_equ ||
             sizemax < 0. ) )
        {
          // copy selected scattering element to temp arrays
          scat_data_tmp.push_back(scat_data[i_ss][i_se]);
          scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
        }
    }
  else if ( sizeparam=="diameter_area_equ_aerodynamical" )
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
    {
      if ( scat_meta[i_ss][i_se].diameter_area_equ_aerodynamical
             > sizemin-sizemin*tolerance &&
           ( sizemax+sizemax*tolerance
               > scat_meta[i_ss][i_se].diameter_area_equ_aerodynamical ||
             sizemax < 0. ) )
        {
          // copy selected scattering element to temp arrays
          scat_data_tmp.push_back(scat_data[i_ss][i_se]);
          scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
        }
    }
  else
    {
      ostringstream os;
      os << "Size parameter " << sizeparam << "is unknown.";
      throw runtime_error ( os.str() );
    }

  // To use a particle species field without associated scattering element
  // data poses a high risk of accidentially neglecting these species. That's
  // unlikely what the user intends. Hence throw error.
  if (scat_meta_tmp.nelem()<1)
    {
      ostringstream os;
      os << "For scattering species " << species << " no scattering "
         << "element matching the requested size range found.\n"
         << "Check scat_data and scat_meta input as well as your size limit "
         << "selection!";
      throw runtime_error ( os.str() );
      }

  scat_meta[i_ss] = scat_meta_tmp;
  scat_data[i_ss] = scat_data_tmp;

  // check if array is empty. should never apply (since we checked the re-worked
  // data before and that error should also catch cases that are empty from the
  // beginning).
  assert ( TotalNumberOfElements(scat_meta) );

}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesExtendTemperature( //WS Output:
                     ArrayOfArrayOfSingleScatteringData& scat_data,
                     // Keywords:
                     const ArrayOfString& scat_species,
                     const String& species,
                     const String& scat_species_delim,
                     const Numeric& T_low,
                     const Numeric& T_high,
                     const Verbosity& )
{
  const bool do_tl = ( T_low >= 0. );
  const bool do_th = ( T_high >= 0. );

  if( do_tl || do_th )
  {
    Index i_ss=-1;
    if( species=="" )
      {
        i_ss = scat_data.nelem()-1;
        if (i_ss==-1)
          {
            ostringstream os;
            os << "No scat_data available. Can not extend temperature range on "
               << "inexistent data.";
            throw runtime_error ( os.str() );
          }
      }
    else
      {
        // first check that sizes of scat_species and scat_data agree
        Index nspecies = scat_species.nelem();
        if ( nspecies != scat_data.nelem() )
          {
            ostringstream os;
            os << "Number of scattering species specified by scat_species does\n"
               << "not agree with number of scattering species in scat_data:\n"
               << "scat_species has " << nspecies << " entries, while scat_data has "
               << scat_data.nelem() << ".";
            throw runtime_error ( os.str() );
          }
        String partfield_name;
        //find the species to handle: compare 'species' to 'partfield' part of
        //scat_species tags
        for ( Index i=0; i<scat_species.nelem(); i++ )
          {
            parse_partfield_name( partfield_name, scat_species[i], scat_species_delim);
            if ( partfield_name==species )
              i_ss=i;
          }
        if ( i_ss<0 )
          {
            ostringstream os;
            os << "Scattering species " << species << " not found among scat_species.";
            throw runtime_error ( os.str() );
          }
      }

    for ( Index i_se=0; i_se<scat_data[i_ss].nelem(); i_se++ )
      {
        const SingleScatteringData ssdo = scat_data[i_ss][i_se];
        const Index nTo = ssdo.T_grid.nelem();
        Index nTn = nTo;
        bool do_htl, do_hth;
        if( nTo > 1 )
          {
            do_htl = ( do_tl && (T_low < ssdo.T_grid[0]) );
            do_hth = ( do_th && (T_high > last(ssdo.T_grid)) );
          }
        else
          {
            do_htl = false;
            do_hth = false;
          }

        if( do_htl || do_hth )
        {
          // create new instance of SingleScatteringData
          SingleScatteringData ssdn;
          Index iToff;

          // determine new temperature grid
          if( do_htl )
            nTn += 1;
          if( do_hth )
            nTn += 1;
          Vector T_grid_new(nTn);
          if( do_htl )
            {
              T_grid_new[0] = T_low;
              iToff = 1;
            }
          else
            {
              iToff = 0;
            }
          for( Index iT=0; iT<nTo; iT++ )
              T_grid_new[iT+iToff] = scat_data[i_ss][i_se].T_grid[iT];
          if( do_hth )
              T_grid_new[nTo+iToff] = T_high;
          ssdn.T_grid = T_grid_new;

          // copy grids and other descriptive data that is to remain identical
          ssdn.ptype = ssdo.ptype;
          ostringstream description;
          description << ssdo.description; // here just copy. we append further
                                           // info below if applicable.
          ssdn.f_grid = ssdo.f_grid;
          ssdn.za_grid = ssdo.za_grid;
          ssdn.aa_grid = ssdo.aa_grid;

          // determine size of current optical property data
          const Index nf = ssdo.f_grid.nelem();
          const Index nzas = ssdo.pha_mat_data.nshelves();
          const Index naas = ssdo.pha_mat_data.nbooks();
          const Index nzai = ssdo.pha_mat_data.npages();
          const Index naai = ssdo.pha_mat_data.nrows();
          const Index nmep = ssdo.pha_mat_data.ncols();
          const Index nmee = ssdo.ext_mat_data.ncols();
          const Index nvea = ssdo.abs_vec_data.ncols();

          // create containers for extended optical property data
          ssdn.pha_mat_data.resize(nf, nTn, nzas, naas, nzai, naai, nmep);
          ssdn.ext_mat_data.resize(nf, nTn, nzai, naai, nmee);
          ssdn.abs_vec_data.resize(nf, nTn, nzai, naai, nvea);

          // copy optical property data
          for( Index iT=0; iT<nTo; iT++ )
            {
              ssdn.pha_mat_data(joker,iT+iToff,joker,joker,joker,joker,joker)
                = ssdo.pha_mat_data(joker,iT,joker,joker,joker,joker,joker);
              ssdn.ext_mat_data(joker,iT+iToff,joker,joker,joker)
                = ssdo.ext_mat_data(joker,iT,joker,joker,joker);
              ssdn.abs_vec_data(joker,iT+iToff,joker,joker,joker)
                = ssdo.abs_vec_data(joker,iT,joker,joker,joker);
            }

          // duplicate optical property data on T-edges if applicable
          if( do_htl )
            {
              ssdn.pha_mat_data(joker,0,joker,joker,joker,joker,joker)
                = ssdn.pha_mat_data(joker,1,joker,joker,joker,joker,joker);
              ssdn.ext_mat_data(joker,0,joker,joker,joker)
                = ssdn.ext_mat_data(joker,1,joker,joker,joker);
              ssdn.abs_vec_data(joker,0,joker,joker,joker)
                = ssdn.abs_vec_data(joker,1,joker,joker,joker);
              description << "\n" << "Low temperature limit extended by"
                          << " duplicating previous low temperature limit"
                          << " single scattering properties.";
            }
          if( do_hth )
            {
              ssdn.pha_mat_data(joker,nTn-1,joker,joker,joker,joker,joker)
                = ssdn.pha_mat_data(joker,nTn-2,joker,joker,joker,joker,joker);
              ssdn.ext_mat_data(joker,nTn-1,joker,joker,joker)
                = ssdn.ext_mat_data(joker,nTn-2,joker,joker,joker);
              ssdn.abs_vec_data(joker,nTn-1,joker,joker,joker)
                = ssdn.abs_vec_data(joker,nTn-2,joker,joker,joker);
              description << "\n" << "High temperature limit extended by"
                          << " duplicating previous high temperature limit"
                          << " single scattering properties.";
            }
          ssdn.description = description.str();
          scat_data[i_ss][i_se] = ssdn;
        }
      }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaDataSingleCategory(
        Matrix&                    particle_masses,
  const ArrayOfArrayOfScatteringMetaData& scat_meta,
  const Verbosity&)
{
  const Index np_total = TotalNumberOfElements(scat_meta);

  particle_masses.resize(np_total,1);

  Index i_se_flat = 0;
  for( Index i_ss=0; i_ss<scat_meta.nelem(); i_ss++ )
  {
      for( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
      {
          if( isnan(scat_meta[i_ss][i_se].mass) ||
              scat_meta[i_ss][i_se].mass <= 0 ||
              scat_meta[i_ss][i_se].mass > 1. )
          {
              ostringstream os;
              os << "A presumably incorrect value found for "
              << "scat_meta[" << i_ss << "][" << i_se << "].mass.\n"
              << "The value is " << scat_meta[i_ss][i_se].mass;
              throw std::runtime_error(os.str());
          }

          if( scat_meta[i_ss][i_se].diameter_volume_equ <= 0 ||
             scat_meta[i_ss][i_se].diameter_volume_equ > 0.5 )
          {
              ostringstream os;
              os << "A presumably incorrect value found for "
              << "scat_meta[" << i_ss << "][" << i_se << "].diameter_volume_equ.\n"
              << "The value is " << scat_meta[i_ss][i_se].diameter_volume_equ;
              throw std::runtime_error(os.str());
          }

          particle_masses(i_se_flat,0) = scat_meta[i_ss][i_se].diameter_volume_equ;

          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaData
                        (//WS Output:
                         Matrix& particle_masses,
                         // WS Input:
                         const ArrayOfArrayOfScatteringMetaData& scat_meta,
                         const Verbosity& )
{
  // resize particle_masses to required diemsions and properly initialize values
  particle_masses.resize ( TotalNumberOfElements(scat_meta), scat_meta.nelem() );
  particle_masses = 0.;

  // calculate and set particle_masses
  Index i_se_flat = 0;
  for ( Index i_ss=0; i_ss<scat_meta.nelem(); i_ss++ )
  {
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
    {
      particle_masses (i_se_flat, i_ss) = scat_meta[i_ss][i_se].diameter_volume_equ;
      i_se_flat++;
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFrompnd_field_raw(//WS Output:
                   Tensor4& pnd_field,
                   //WS Input
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const ArrayOfGriddedField3& pnd_field_raw,
                   const Index& atmosphere_dim,
                   const ArrayOfIndex& cloudbox_limits,
                   const Index& zeropadding,
                   const Verbosity& verbosity)
{
  // Basic checks of input variables
  //
  // Particle number density data
  // 
  if (pnd_field_raw.empty())
  {
    ostringstream os;
    os << "No particle number density data given. Please use WSMs\n"
       << "*ParticleTypeInit* and *ParticleTypeAdd(All)* for reading\n"
       << "scattering element data.\n";
    throw runtime_error(os.str());
  }
  
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  ArrayOfIndex cloudbox_limits_tmp;
  /*if ( cloudbox_limits.empty() )
    {
      //If no limits set, the cloud(box) is supposed to cover the
      //complete atmosphere. This particularly to facilitate use of
      //scat_data_single&pnd_fields for non-scatt, greybody cloud calculations.
      //Hence, set the limits to first & last elements of the different grids.
      //Note: no margin left at lat/lon_grid edges!.
      cloudbox_limits_tmp.resize(2*atmosphere_dim);

      // Pressure limits
      cloudbox_limits_tmp[0] = 0;
      cloudbox_limits_tmp[1] = p_grid.nelem() - 1;
      // Latitude limits
      if( atmosphere_dim >= 2 )
        {
          cloudbox_limits_tmp[2] = 0;
          cloudbox_limits_tmp[3] = lat_grid.nelem() - 1;
        }
      // Longitude limits
      if( atmosphere_dim == 3 )
        {
          cloudbox_limits_tmp[4] = 0;
          cloudbox_limits_tmp[5] = lon_grid.nelem() - 1;
        }
    }
  else */
  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");
  else
    cloudbox_limits_tmp = cloudbox_limits;

  // Check that pnd_field_raw has at least 2 grid-points in each dimension.
  // Otherwise, interpolation further down will fail with assertion.
  
  for (Index d = 0; d < atmosphere_dim; d++)
    {
      for (Index i = 0; i < pnd_field_raw.nelem(); i++)
        {
          if (pnd_field_raw[i].get_grid_size(d) < 2)
            {
              ostringstream os;
              os << "Error in pnd_field_raw data. ";
              os << "Dimension " << d << " (name: \"";
              os << pnd_field_raw[i].get_grid_name(d);
              os << "\") has only ";
              os << pnd_field_raw[i].get_grid_size(d);
              os << " element";
              os << ((pnd_field_raw[i].get_grid_size(d)==1) ? "" : "s");
              os << ". Must be at least 2.";
              throw runtime_error(os.str());
            }
        }
    }
  const Index Np_cloud = cloudbox_limits_tmp[1]-cloudbox_limits_tmp[0]+1;
  
  ConstVectorView p_grid_cloud = p_grid[Range(cloudbox_limits_tmp[0], Np_cloud)];

  // Check that no scatterers exist outside the cloudbox
  chk_pnd_field_raw_only_in_cloudbox(atmosphere_dim, pnd_field_raw,
                                     p_grid, lat_grid, lon_grid,
                                     cloudbox_limits_tmp);

  //==========================================================================
  if ( atmosphere_dim == 1)
    {
        ArrayOfGriddedField3 pnd_field_tmp;

        GriddedFieldPRegrid(pnd_field_tmp, p_grid_cloud, pnd_field_raw,
                            1, zeropadding, verbosity);

        FieldFromGriddedField(pnd_field,
                              p_grid_cloud,
                              pnd_field_tmp[0].get_numeric_grid(1),
                              pnd_field_tmp[0].get_numeric_grid(2),
                              pnd_field_tmp, verbosity);
    }
  else if(atmosphere_dim == 2)
    {
      const Index Nlat_cloud = cloudbox_limits_tmp[3]-cloudbox_limits_tmp[2]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits_tmp[2],Nlat_cloud)];           

      if (zeropadding)
        {
          // FIXME: OLE: Implement this
          CREATE_OUT0;
          out0 << "WARNING: zeropadding currently only supported for 1D.";
        }

      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 1 );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the scattering element:
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
                  pnd_field_raw[i].data(joker,joker,0), gp_p, gp_lat );
        }
    }
  else
    {
      const Index Nlat_cloud = cloudbox_limits_tmp[3]-cloudbox_limits_tmp[2]+1;
      const Index Nlon_cloud = cloudbox_limits_tmp[5]-cloudbox_limits_tmp[4]+1;

      if (zeropadding)
        {
          // FIXME: OLE: Implement this
          CREATE_OUT0;
          out0 << "WARNING: zeropadding currently only supported for 1D.";
        }

      ConstVectorView lat_grid_cloud =
        lat_grid[Range(cloudbox_limits_tmp[2],Nlat_cloud)];           
      ConstVectorView lon_grid_cloud = 
        lon_grid[Range(cloudbox_limits_tmp[4],Nlon_cloud)];
      
      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 
                        Nlon_cloud );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      ArrayOfGridPos gp_lon(Nlon_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the scattering element types:
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
                  pnd_field_raw[i].data, gp_p, gp_lat, gp_lon );
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldExpand1D(Tensor4&        pnd_field,
                       const Index&    atmosphere_dim,
                       const Index&    cloudbox_on,    
                       const ArrayOfIndex&   cloudbox_limits,
                       const Index&    nzero,
                       const Verbosity&)
{
/*  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
*/
  if( atmosphere_dim == 1 )
    { throw runtime_error( "No use in calling this method for 1D." ); }
  if( !cloudbox_on )
    { throw runtime_error( 
                "No use in calling this method with cloudbox off." ); }

  if( nzero < 1 )
    { throw runtime_error( "The argument *nzero* must be > 0." ); }

  // Sizes
  const Index   npart = pnd_field.nbooks();
  const Index   np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index   nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
        Index   nlon = 1;
  if( atmosphere_dim == 3 )  
    { nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1; }

  if( pnd_field.npages() != np  ||  pnd_field.nrows() != 1  ||  
      pnd_field.ncols() != 1 )
    { throw runtime_error( "The input *pnd_field* is either not 1D or does not "
                           "match pressure size of cloudbox." );}

  // Temporary container
  Tensor4 pnd_temp = pnd_field;

  // Resize and fill
  pnd_field.resize( npart, np, nlat, nlon );
  pnd_field = 0;
  //
  for( Index ilon=nzero; ilon<nlon-nzero; ilon++ )
    {
      for( Index ilat=nzero; ilat<nlat-nzero; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for( Index is=0; is<npart; is++ )
                { pnd_field(is,ip,ilat,ilon) = pnd_temp(is,ip,0,0); }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldZero(//WS Output:
                   Tensor4& pnd_field,
                   ArrayOfArrayOfSingleScatteringData& scat_data,
                   //WS Input:
                   const Index& atmosphere_dim,
                   const ArrayOfIndex& cloudbox_limits,
                   const Verbosity&)
{
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  
  //Resize pnd_field and set it to 0:
  Index np = cloudbox_limits[1]-cloudbox_limits[0]+1;
  Index nlat=1, nlon=1;
  if( atmosphere_dim > 1 )
    {
      nlat = cloudbox_limits[3]-cloudbox_limits[2]+1;
      if( atmosphere_dim > 2 )
        {
          nlon = cloudbox_limits[5]-cloudbox_limits[4]+1;
        }
    }

    
  // Do only reset scat_data if it has not been set yet.
  // There's no need otherwise, and it's rather unpractical for testing when
  // doing so: we might want to do some actual calcs with the scat_data later
  // on. So why throw it away?
  const Index N_se = TotalNumberOfElements(scat_data);
  if( N_se > 0 )
    {
      pnd_field.resize(N_se, np, nlat, nlon);
    }
  else
    {
      pnd_field.resize(1, np, nlat, nlon);
     
      //Resize scat_data and set it to 0:
      // Number of scattering elements
      scat_data.resize(1);
      scat_data[0].resize(1);
      scat_data[0][0].ptype = PTYPE_TOTAL_RND;
      scat_data[0][0].description = " ";
      // Grids which contain full ranges which one wants to calculate
      nlinspace(scat_data[0][0].f_grid, 1e9, 3.848043e+13 , 5);
      nlinspace(scat_data[0][0].T_grid, 0, 400, 5);
      nlinspace(scat_data[0][0].za_grid, 0, 180, 5);
      nlinspace(scat_data[0][0].aa_grid, 0, 360, 5);
      // Resize the data arrays
      scat_data[0][0].pha_mat_data.resize(5,5,5,1,1,1,6);
      scat_data[0][0].pha_mat_data = 0.;
      scat_data[0][0].ext_mat_data.resize(5,5,1,1,1);
      scat_data[0][0].ext_mat_data = 0.;
      scat_data[0][0].abs_vec_data.resize(5,5,1,1,1);
      scat_data[0][0].abs_vec_data = 0.;
    }

  pnd_field = 0.;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFromscat_speciesFields (//WS Output:
                     Tensor4& pnd_field,
                     //WS Input:
                     const Index& atmosphere_dim,
                     const Index& cloudbox_on,
                     const ArrayOfIndex& cloudbox_limits,
                     const Tensor4& scat_species_mass_density_field,
                     const Tensor4& scat_species_mass_flux_field,
                     const Tensor4& scat_species_number_density_field,
                     const Tensor4& scat_species_mean_mass_field _U_,
                     const Tensor3& t_field,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const ArrayOfString& scat_species,
                     const String& delim,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;

  // Cloudbox on/off?
  if ( !cloudbox_on )
  {
    /* Must initialise pnd_field anyway; but empty */
    pnd_field.resize(0, 0, 0, 0);

    CREATE_OUT0;
    out0 << "  Cloudbox is off, pnd_field is set to empty.\n";
    return;
  }

  // Check consistency of scat_species and scat_data
  if ( scat_species.nelem() != scat_meta.nelem() )
    {
      ostringstream os;
      os << "Number of scattering species specified by scat_species does\n"
         << "not agree with number of scattering species in scat_meta:\n"
         << "scat_species has " << scat_species.nelem()
         << " entries, while scat_meta has " << scat_meta.nelem() << ".";
      throw runtime_error ( os.str() );
    }

  // ------- set pnd_field boundaries to cloudbox boundaries -------------------
  //set pnd_field boundaries
  ArrayOfIndex limits(6);
  //pressure
  limits[0] = cloudbox_limits[0];
  limits[1] = cloudbox_limits[1]+1;
  //latitude
  if ( atmosphere_dim > 1 )
  {
      limits[2] = cloudbox_limits[2];
      limits[3] = cloudbox_limits[3]+1;
  }
  else
  {
      limits[2] = 0;
      limits[3] = 1;
  }
  if ( atmosphere_dim > 2 )
  {
      limits[4] = cloudbox_limits[4];
      limits[5] = cloudbox_limits[5]+1;
  }
  else
  {
      limits[4] = 0;
      limits[5] = 1;
  }

  /* Do some checks. Not foolproof, but catches at least some. */

  if ((limits[1] > scat_species_mass_density_field.npages()) ||
      (limits[1] > t_field.npages()) ||
      (limits[3] > scat_species_mass_density_field.nrows()) ||
      (limits[3] > t_field.nrows()) ||
      (limits[5] > scat_species_mass_density_field.ncols()) ||
      (limits[5] > t_field.ncols()))
  {
    ostringstream os;
    os << "Cloudbox out of bounds compared to fields. "
       << "Upper limits: (p, lat, lon): "
       << "(" << limits[1] << ", " << limits[3] << ", " << limits[5] << "). "
       << "*scat_species_mass_density_field*: "
       << "(" << scat_species_mass_density_field.npages() << ", "
       << scat_species_mass_density_field.nrows() << ", "
       << scat_species_mass_density_field.ncols() << "). "
       << "*t_field*: "
       << "(" << t_field.npages() << ", "
       << t_field.nrows() << ", "
       << t_field.ncols() << ").";
    throw runtime_error(os.str());
  }

  //resize pnd_field to required atmospheric dimension and scattering elements
  pnd_field.resize ( TotalNumberOfElements(scat_meta), limits[1]-limits[0],
                     limits[3]-limits[2], limits[5]-limits[4] );
  ArrayOfIndex intarr;

  //-------- Start pnd_field calculations---------------------------------------

  // loop over nelem of scat_species
  for ( Index i_ss=0; i_ss<scat_species.nelem(); i_ss++ )
  {

    //cout << "for scat species #" << i_ss << " using PSD ";
    String psd_param;
    String partfield_name;
    String psd;

    //split String and copy to ArrayOfString
    parse_psd_param( psd_param, scat_species[i_ss], delim);
    parse_partfield_name( partfield_name, scat_species[i_ss], delim);

    if (scat_meta[i_ss].nelem() == 0)
    {
        ostringstream os;
        os << "*scat_meta* for scattering species '" << scat_species[i_ss] << "' is empty.";
        throw std::runtime_error(os.str());
    }

  /*
   * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * NOTE: when adding further distributions here, document them (particularly
   * the tags) in the WSM online documentation in methods.cc. Also, create a
   * wrapper WSM to the dNdD core method.
   * See ARTS Developer Guide for details.
   * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   */

    if ( psd_param == "mergedpsd" )
    {
      ostringstream os;
      os << "Scattering species #" << i_ss << " data seems to originate from "
         << "*ScatSpeciesMerge*,\n"
         << "This has unrevertibly modified *scat_data*, *scat_meta*, and "
         << "*scat_species*,\n"
         << "and destroyed their links to the scat_species_*_field data.";
      throw runtime_error ( os.str() );
    }

    //---- pnd_field calculations for MH97 -------------------------------
    if ( psd_param.substr(0,4) == "MH97" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should should only be applied to cloud"
                 << " ice.\n";
        }

        pnd_fieldMH97( pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }

    //---- pnd_field calculations for H11 ----------------------------
    else if ( psd_param == "H11" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        pnd_fieldH11 (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }
    
    //---- pnd_field calculations for H13 ----------------------------
    else if ( psd_param == "H13" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC"  )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        pnd_fieldH13 (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }
    
//---- pnd_field calculations for H13Shape ----------------------------
    else if ( psd_param == "H13Shape" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC"  )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        pnd_fieldH13Shape (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }
      
    //---- pnd_field calculations for F07TR ----------------------------
    else if ( psd_param == "F07TR" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to cloud or precipitating"
            << " ice.\n";
        }
        
        pnd_fieldF07 (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      t_field, "TR", limits,
                      scat_meta, i_ss, scat_species[i_ss], delim, verbosity);
    }
      
      //---- pnd_field calculations for F07ML ----------------------------
    else if ( psd_param == "F07ML" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to cloud or precipitating"
            << " ice.\n";
        }
        
        pnd_fieldF07 (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      t_field, "ML", limits,
                      scat_meta, i_ss, scat_species[i_ss], delim, verbosity);
    }
      
      //---- pnd_field calculations for S2M_LWC ----------------------------
    else if ( psd_param == "S2M_LWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to liquid cloud water.\n";
        }
        
        pnd_fieldS2M (pnd_field,
                        scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                        scat_species_number_density_field ( i_ss, joker, joker, joker ),
                        limits,
                        scat_meta,
                        i_ss,
                        scat_species[i_ss],
                        delim,
                        verbosity);
    }
      
      //---- pnd_field calculations for S2M_LWC_M ----------------------------
    else if ( psd_param == "S2M_LWC_M" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
            << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to liquid cloud water.\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_mean_mass_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_IWC ----------------------------
    else if ( psd_param == "S2M_IWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param<< ".\n"
            << psd_param << " should only be applied to cloud ice.\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      //---- pnd_field calculations for S2M_IWC_M ----------------------------
    else if ( psd_param == "S2M_IWC_M" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
            << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to cloud ice.\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_mean_mass_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_RWC ----------------------------
    else if ( psd_param == "S2M_RWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Rain" && partfield_name != "RWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to precipitating liquid water\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_RWC_M ----------------------------
    else if ( psd_param == "S2M_RWC_M" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Rain" && partfield_name != "RWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
            << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to precipitating liquid water\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_mean_mass_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_SWC ----------------------------
    else if ( psd_param == "S2M_SWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to snow\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }

      //---- pnd_field calculations for S2M_SWC_M ----------------------------
    else if ( psd_param == "S2M_SWC_M" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
            << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to snow\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_mean_mass_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_GWC ----------------------------
    else if ( psd_param == "S2M_GWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Graupel" && partfield_name != "GWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to graupel\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_GWC_M ----------------------------
    else if ( psd_param == "S2M_GWC_M" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Graupel" && partfield_name != "GWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
            << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to graupel\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_mean_mass_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_HWC ----------------------------
    else if ( psd_param == "S2M_HWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Hail" && partfield_name != "HWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to hail\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_HWC_M ----------------------------
    else if ( psd_param == "S2M_HWC_M" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "Hail" && partfield_name != "HWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
            << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to hail\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_mean_mass_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for MGD_LWC ----------------------------
    else if ( psd_param == "MGD_LWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to liquid cloud water.\n";
        }
        
        pnd_fieldMGD_LWC (pnd_field,
                        scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                        limits,
                        scat_meta, i_ss,
                        scat_species[i_ss], delim,
                        verbosity);
    }
      
      //---- pnd_field calculations for MGD_LWC ----------------------------
    else if ( psd_param == "MGD_IWC" )
    {
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
            << psd_param << " should only be applied to cloud ice.\n";
        }
        
        pnd_fieldMGD_IWC (pnd_field,
                          scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                          limits,
                          scat_meta, i_ss,
                          scat_species[i_ss], delim,
                          verbosity);
    }
      
    //---- pnd_field calculations for W16 -------------------------------
    else if ( psd_param == "W16" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "Rain" && partfield_name != "RWC" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should only be applied to liquid"
                 << " precipitation (aka rain).\n";
        }

        pnd_fieldW16( pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta, i_ss,
                      scat_species[i_ss], delim,
                      verbosity);
    }

    //---- pnd_field calculations for MP48 -------------------------------
    else if ( psd_param == "MP48" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "Rain" && partfield_name != "Snow"
            && partfield_name != "RR" && partfield_name != "SR" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should only be applied to liquid or frozen"
                 << " precipitation.\n";
        }

        pnd_fieldMP48( pnd_field,
                       scat_species_mass_flux_field ( i_ss, joker, joker, joker ),
                       limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }

    // ---- pnd_field calculations for liquid water clouds ---------------------
    //      (parameters from Hess98, Stratus continental)
    else if ( psd_param == "H98_STCO" || psd_param == "STCO" || psd_param == "liquid" )
    {

        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd_param << ".\n"
                 << psd_param << " should should only be applied to liquid or frozen"
                 << " precipitation.\n";
        }

        pnd_fieldH98 (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }

    else
    {
        ostringstream os;
        os << "Size distribution function '" << psd_param << "' is unknown!\n";
        throw runtime_error( os.str() );
    }

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * NOTE: when adding further distributions here, document them (particularly
   * the tags) in the WSM online documentation in methods.cc. Also, create a
   * wrapper WSM to the dNdD core method.
   * See ARTS Developer Guide for details.
   * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromdNdD (//WS Output:
              Vector& pnd,
              //WS Input:
              const Vector& dNdD,
              const Vector& diameter,
              const Numeric& total_content,
              const Vector& scatelem_content,
              const Verbosity& verbosity)
{
  pnd.resize( diameter.nelem() );

  // scale dNdD by size bin width
  if (diameter.nelem() > 1)
      bin_integral( pnd, diameter, dNdD ); //[# m^-3]
  else
      pnd = dNdD;

  // scaling pnd to real mass/number density/flux (some PSDs have implicit
  // scaling - then this is only a check -, others don't
  chk_pndsum ( pnd, total_content, scatelem_content,
               0, 0, 0, "your field", verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_size_gridFromScatMeta(
          Vector&                             pnd_size_grid,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&                              scat_index,
    const String&                             unit,
    const Verbosity& )
{
  // Sizes
  const Index nss = scat_meta.nelem();

  // Checks
  if( nss == 0 )
    throw runtime_error( "*scat_meta* is empty!" );
  if( nss < scat_index+1 )
    {
      ostringstream os;
      os << "Selected scattering species index is " << scat_index << " but this "
         << "is not allowed since *scat_meta* has only" << scat_meta.nelem()
         << " elements.";
      throw runtime_error(os.str());
    }

  // Create size grid
  //
  const Index nse = scat_meta[scat_index].nelem();
  //
  pnd_size_grid.resize( nse );
  //
  if( unit == "dveq" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_volume_equ; }
    }
  else if( unit == "dmax" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_max; }
    }
  else if( unit == "mass" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].mass; }
    }
  else if( unit == "area" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_area_equ_aerodynamical; }
    }
  else
    {
      ostringstream os;
      os << "You have selected the unit: " << unit 
         << "while accepted choices are: \"dveq\", \"dmax\", \"mass\" and \"area\"";
      throw runtime_error(os.str());
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsdBasic(
         Matrix&    pnd_data,
         Tensor3&   dpnd_data_dx,
   const Vector&    pnd_size_grid,
   const Matrix&    psd_data,
   const Vector&    psd_size_grid,
   const Tensor3&   dpsd_data_dx,
   const Verbosity& )
{
  // This code is just for testing. Shall be replaced with code using
  // piece-wise linear representation when mapping psd to pnd.
  
  // Some sizes 
  const Index np  = psd_data.nrows();
  const Index ng  = psd_size_grid.nelem();
        Index ndx = 0;
  const bool  do_dx = !dpsd_data_dx.empty();

  // Checks
  if( ng < 2 )
    throw runtime_error( "The method requires that length of *psd_size_grid* is >= 2." );
  if( ng != pnd_size_grid.nelem() )
    throw runtime_error( "The method requires that *psd_size_grid* and "
                         "*pnd_size_grid* have same length." );
  for( Index i=0; i<ng; i++ )
    {
      if( psd_size_grid[i] != pnd_size_grid[i] )
        throw runtime_error( "So far, the method requires that *psd_size_grid* and "
                             "*pnd_size_grid* are identical." );
    }
  if( psd_data.ncols() != ng )
    throw runtime_error( "Number of columns in *psd_data* and length of "
                         "*psd_size_grid* must match." );
  if( do_dx )
    {
      if( dpsd_data_dx.ncols() != ng )
        throw runtime_error( "Number of columns in *dpsd_data_dx* and length of "
                             "*psd_size_grid* must match." );
      ndx = dpsd_data_dx.npages();
      dpnd_data_dx.resize( ndx, np, ng );
    }
  else
    { dpnd_data_dx.resize( 0, 0, 0 ); }

  if( !is_increasing( psd_size_grid ) )
    throw runtime_error( "*psd_size_grid* must be strictly increasing." );    


  // Also, is my scheme to define bin sizes the same as you used?
  // Asking as I get deviating results, and I suspect it could come from this
  // part 

  // dpnd_dx is sized above    
  pnd_data.resize( np, ng );

  
  // Calculate
  Numeric binsize;
  for ( Index i=0; i<ng; i++ )
    {
      // This bin is twice the half-distance to point 1, but could be limited
      // by 0 in lower end
      if( i == 0 )
        {
          const Numeric dd = ( psd_size_grid[1] - psd_size_grid[0] ) / 2.0;
          binsize = dd + min( dd, psd_size_grid[0] );
        }
      // This bin is twice the half-distance to closest point      
      else if( i == ng-1 )
        { binsize = 2.0 * ( psd_size_grid[i] - psd_size_grid[i-1] ) / 2.0; }
      // This bin is the sum of the two half-distances      
      else
        { binsize = ( psd_size_grid[i+1] - psd_size_grid[i-1] ) / 2.0; }

      for( Index ip=0; ip<np; ip++ )
        { pnd_data(ip,i) = binsize * psd_data(ip,i); }

      if( do_dx )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for ( Index ix=0; ix<ndx; ix++ )
                { dpnd_data_dx(ix,ip,i) = binsize * dpsd_data_dx(ix,ip,i); }
            }
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsdBasic2(
         Matrix&    pnd_data,
         Tensor3&   dpnd_data_dx,
   const Vector&    pnd_size_grid,
   const Matrix&    psd_data,
   const Vector&    psd_size_grid,
   const Tensor3&   dpsd_data_dx,
   const Index&     quad_order,
   const Verbosity& )
{
  // Some sizes 
  const Index np  = psd_data.nrows();
  const Index ng  = psd_size_grid.nelem();
        Index ndx = 0;
  const bool  do_dx = !dpsd_data_dx.empty();

  // Checks
  if( ng < 2 )
    throw runtime_error( "The method requires that length of *psd_size_grid* is >= 2." );
  if( ng != pnd_size_grid.nelem() )
    throw runtime_error( "So far, the method requires that *psd_size_grid* and "
                         "*pnd_size_grid* have same length." );
  for( Index i=0; i<ng; i++ )
    {
      if( psd_size_grid[i] != pnd_size_grid[i] )
        throw runtime_error( "So far, the method requires that *psd_size_grid* and "
                             "*pnd_size_grid* are identical." );
    }
  if( psd_data.ncols() != ng )
    throw runtime_error( "Number of columns in *psd_data* and length of "
                         "*psd_size_grid* must match." );
  if( do_dx )
    {
      if( dpsd_data_dx.ncols() != ng )
        throw runtime_error( "Number of columns in *dpsd_data_dx* and length of "
                             "*psd_size_grid* must match." );
      ndx = dpsd_data_dx.npages();
      dpnd_data_dx.resize( ndx, np, ng );
    }
  else
    { dpnd_data_dx.resize( 0, 0, 0 ); }

  // the psd_data calculation methods allow unsorted size grids. so, we should
  // either provide a WSM that sorts size grids or do sorting here (note,
  // particularly there is no guarantee that sizes in scat_meta are sorted. a
  // special problem there is, that sorting might depend on the size parameter
  // chosen.
  if( !is_increasing( psd_size_grid ) )
    throw runtime_error( "*psd_size_grid* must be strictly increasing." );    

  // dpnd_dx is sized above    
  pnd_data.resize( np, ng );

  // Calculate
  Vector quadweights( ng );
  bin_quadweights( quadweights, psd_size_grid, quad_order );

  for ( Index i=0; i<ng; i++ )
    {
      for( Index ip=0; ip<np; ip++ )
        { pnd_data(ip,i) = quadweights[i] * psd_data(ip,i); }

      if( do_dx )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for ( Index ix=0; ix<ndx; ix++ )
                { dpnd_data_dx(ix,ip,i) = quadweights[i] *
                                          dpsd_data_dx(ix,ip,i); }
            }
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pndAdjustFromScatMeta(
          Matrix&                             pnd_data,
//          Tensor3&                            dpnd_data_dx,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const String&                             pnd_agenda_input_tag,
    const Index&                              scat_index,
    const Verbosity&                          verbosity )
{
  const Index np  = pnd_data.nrows();
  if( pnd_agenda_input.nrows() != np )
    throw runtime_error( "Number of rows in *pnd_data* and *pnd_agenda_input*"
                         " must be equal." );

  const Index nin = pnd_agenda_input_names.nelem();
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );

  const Index nse  = pnd_data.ncols();
  if( scat_meta[scat_index].nelem() != nse )
    {
      ostringstream os;
      os << "Number of columns in *pnd_data* and elements in *scat_meta*\n"
         << "of species " << scat_index << " must be equal.";
      throw runtime_error(os.str());
    }

  Index not_found = 1;
  Index fieldID = 0;
  while( not_found && (fieldID<nin) )
    {
      if( pnd_agenda_input_names[fieldID]==pnd_agenda_input_tag )
        not_found = 0;
      else
        fieldID += 1;
    }
  if( not_found )
    {
      ostringstream os;
      os << "Requested adjustment field \"" << pnd_agenda_input_tag
         << "\" not found in *pnd_agenda_input_names*.";
      throw runtime_error(os.str());
    }

  Vector mass( nse );
  // we have to assume correspondence of scattering element order in pnd_data
  // and scat_meta[scat_index]. (any way to ensure that?)
  for ( Index i=0; i<nse; i++ )
    {
      if ( isnan(scat_meta[scat_index][i].mass) )
        {
          ostringstream os;
          os << "No mass data available for scattering element #"
             << i << " of scattering species with index "
             << scat_index << ".";
          throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_index][i].mass;
    }

  for( Index ip=0; ip<np; ip++ )
    {
      // scaling pnd to real mass/number density/flux (some PSDs have implicit
      // scaling - then this should only a minor adjustment, rather a check -,
      // but others don't
      Vector pnd = pnd_data(ip,joker);
      chk_pndsum ( pnd, pnd_agenda_input(ip,fieldID), mass,
                   ip, -1, -1, pnd_agenda_input_tag, verbosity );
      pnd_data(ip,joker) = pnd;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MassSizeParamsFromScatMeta(
          Numeric&                            alpha,
          Numeric&                            beta,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&                              scat_index,
    const Numeric&                            beta_default,
    const Verbosity& )
{
  // Sizes
  const Index nss = scat_meta.nelem();

  // Checks
  if( nss == 0 )
    throw runtime_error( "*scat_meta* is empty!" );
  if( nss < scat_index+1 )
    {
      ostringstream os;
      os << "Selected scattering species index is " << scat_index
         << " but *scat_meta* contains only " << scat_meta.nelem()
         << " species.";
      throw runtime_error(os.str());
    }

  // Extract size and mass data from scat_meta
  const Index nse = scat_meta[scat_index].nelem();

  if( nse==0 )
    {
      ostringstream os;
      os << "*scat_meta* of selected scattering species with index "
         << scat_index << " contains has no scattering elements.";
      throw runtime_error(os.str());
    }

  if( nse>1 )
    {
      ArrayOfIndex intarr;
      Vector dmax_unsorted( nse );
      Vector dmax( nse ), log_D( nse );
      Vector mass( nse ), log_m( nse );
      Vector q;

      for ( Index i=0; i<nse; i++ )
      {
        if ( isnan(scat_meta[scat_index][i].diameter_max) )
        {
            ostringstream os;
            os << "No maximum diameter data available for scattering element #"
               << i << " of scattering species with index " << scat_index
               << ".";
            throw runtime_error( os.str() );
        }
        dmax_unsorted[i] = scat_meta[scat_index][i].diameter_max;
      }

      get_sorted_indexes(intarr, dmax_unsorted);
    
      // extract scattering meta data
      for ( Index i=0; i<nse; i++ )
      {
        dmax[i] = dmax_unsorted[intarr[i]]; // [m]
        
        if ( isnan(scat_meta[scat_index][intarr[i]].mass) )
        {
            ostringstream os;
            os << "No mass data available for scattering element #"
               << intarr[i] << " of scattering species with index "
               << scat_index << ".";
            throw runtime_error( os.str() );
        }
        mass[i] = scat_meta[scat_index][intarr[i]].mass; // [kg]
        
        // logarithm of Dmax, needed for estimating mass-dimension-relationship
        log_D[i] = log(dmax[i]);
        
        // logarithm of the mass, even though it is  little weird to have
        // a logarithm of something with a unit...
        log_m[i] = log(mass[i]);
      }
    
      // estimate mass-dimension relationship from meta data by linear regression
      // Assumption of a power law for the mass dimension relationship
      // Approach: log(m) = log(alpha)+beta*log(dmax/D0)
      linreg(q,log_D, log_m);
      alpha = exp(q[0]);
      beta = q[1];
    }
  else
    {
      // for a monodispersion we can't estimate the m-D relation (1 relation, 2
      // unknowns), hence we fix one of the parameters and calculate the other
      // such that we have them consistent. but shouldn't make any difference on
      // the end result whatever we choose here (all ice has to end up with this
      // scattering anyways)
      beta = beta_default;
      if ( isnan(scat_meta[scat_index][0].diameter_max) ||
           isnan(scat_meta[scat_index][0].mass) )
        {
          ostringstream os;
          os << "Maximum diameter or mass data missing for scattering element #"
               << 0 << " of scattering species with index " << scat_index
               << ".";
          throw runtime_error( os.str() );
        }
      alpha = scat_meta[scat_index][0].mass / 
              pow(scat_meta[scat_index][0].diameter_max,beta);
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdF07 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const String&                             regime,
    const Numeric&                            alpha,
    const Numeric&                            beta,
    const Numeric&                            t_min,
    const Numeric&                            t_max,
    const Numeric&                            t_min_psd,
    const Numeric&                            t_max_psd,
    const Numeric&                            beta_min, 
    const Numeric&                            beta_max, 
    //const Numeric&                            density_deviation,
    const Index&                              picky, 
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 2 )
    throw runtime_error( "*pnd_agenda_input* must have two columns (SWC and Temperature)." );
  if( pnd_agenda_input_names[0] != "SWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"SWC\"." );
  if( pnd_agenda_input_names[1] != "Temperature" )
    throw runtime_error( "Element 1 of *pnd_agenda_input_names* must be "
                         "\"Temperature\"." );
  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "SWC" )
        throw runtime_error( "With F07, the only valid option for "
                             "*dpnd_data_dx_names* is: \"SWC\"." );
    }        
  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // SWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  // Some sanity checks. If not picky, we just ignore the insanities.
  if( picky )
    {
      // Sane beta?
      if( beta < beta_min  ||  beta > beta_max )
        {
          ostringstream os;
          os << "Method called with a mass-dimension-relation exponent beta of "
             << beta << ".\n"
             << "This is outside the specified allowed range: ["
             << beta_min << "," << beta_max << "]";
          throw runtime_error(os.str());
        }

      // Do alpha and beta result in sane mass for ice particles (i.e. density
      // shouldn't be > DENSITY_ICE? we check the smallest and largest
      // particles.

      // first, find the smallest and largest particles D_small and D_large (as
      // psd_size grid is not necessarily sorted).
      Vector D_lims{1e3,-1.}; // 1km should be large enough to find
                                         // one that is smaller...
      for( Index iD=0; iD<nsi; iD++ )
        {
          if( D_lims[0] > psd_size_grid[iD] ) D_lims[0]=psd_size_grid[iD];
          if( D_lims[1] < psd_size_grid[iD] ) D_lims[1]=psd_size_grid[iD];
        }

      // now calc resulting densities for solid spheres D_small and D_large
      // diameters
      /* that seems unpractical. for beta<3., small particles will frequently
         exceed ice density quite significantly.
      for( Index iD=0; iD<2; iD++ )
        {
          Numeric rho = 6.*alpha/PI * pow(D_lims[iD],(beta-3.));
          if( rho*(1.+density_deviation) > DENSITY_OF_ICE )
            {
              ostringstream os;
              os << "Method called with a mass-dimension-relation parameters"
                 << " alpha=" << alpha << " and beta=" << beta << "\n"
                 << "results in minimum density of " << rho << ".\n"
                 << "This is unphysical for ice particles and beyond the"
                 << " allowed unphysicality threshold: "
                 << rho*(1.+density_deviation) << ".";
              throw runtime_error(os.str());
            }
        }
      */
    }

  for( Index ip=0; ip<np; ip++ )
    {
      
      // Extract the input variables
      Numeric swc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input(ip,1);

      // No calc needed if swc==0 and no jacobians requested.
      if( (swc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Negative swc?
      Numeric psd_weight = 1.0;
      if( swc < 0 )
        {
          psd_weight = -1.0;
          swc       *= -1.0;
        }
      
      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: ["
                 << t_min << "," << t_max << "]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
      if( t < t_min_psd )
        { t = t_min_psd; }
      else if( t > t_max_psd )
        { t = t_max_psd; }
    
      // Calculate PSD
      Vector psd_1p(nsi);
      if( swc != 0 )
        {
          psd_snow_F07 ( psd_1p, psd_size_grid, swc, t, alpha, beta, regime );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the latter value becomes the perturbation for IWC=0.
          const Numeric dswc = max( 0.001*swc, 1e-7 );
          const Numeric swcp = swc + dswc;
          psd_snow_F07 ( psd_1p, psd_size_grid, swcp, t, alpha, beta, regime );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / dswc; }
        }   
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void psdMH97 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const Numeric&                            t_min, 
    const Numeric&                            t_max, 
    const Numeric&                            t_min_psd, 
    const Numeric&                            t_max_psd, 
    const Index&                              picky, 
    const Index&                              noisy,
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 2 )
    throw runtime_error( "*pnd_agenda_input* must have two columns (IWC and Temperature)." );
  if( pnd_agenda_input_names[0] != "IWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"IWC\"." );
  if( pnd_agenda_input_names[1] != "Temperature" )
    throw runtime_error( "Element 1 of *pnd_agenda_input_names* must be "
                         "\"Temperature\"." );
  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "IWC" )
        throw runtime_error( "With MH97, the only valid option for "
                             "*dpnd_data_dx_names* is: \"IWC\"." );
    }        
  if( noisy  &&   ndx )
    throw runtime_error( "Jacobian calculations and \"noisy\" can not be "
                         "combined." );

  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // IWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  for( Index ip=0; ip<np; ip++ )
    {
      
      // Extract the input variables
      Numeric iwc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input(ip,1);

      // No calc needed if iwc==0 and no jacobians requested.
      if( (iwc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Negative iwc?
      Numeric psd_weight = 1.0;
      if( iwc < 0 )
        {
          psd_weight = -1.0;
          iwc       *= -1.0;
        }
      
      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: ["
                 << t_min << "," << t_max << "]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      if( t < t_min_psd )
        { t = t_min_psd; }
      else if( t > t_max_psd )
        { t = t_max_psd; }
  
      // Calculate PSD
      Vector psd_1p(nsi);
      if( iwc != 0 )
        {
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwc, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the latter value becomes the perturbation for IWC=0.
          const Numeric diwc = max( 0.001*iwc, 1e-7 );
          const Numeric iwcp = iwc + diwc;
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwcp, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / diwc; }
        }   
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void psdW16 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const Index&                              picky _U_, 
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 1 )
    throw runtime_error( "*pnd_agenda_input* must have one column." );
  if( pnd_agenda_input_names[0] != "RWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"RWC\"." );
  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "RWC" )
        throw runtime_error( "With W16, the only valid option for "
                             "*dpnd_data_dx_names* is: \"RWC\"." );
    }        

  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // RWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  for( Index ip=0; ip<np; ip++ )
    {
      
      // Extract the input variables
      Numeric rwc = pnd_agenda_input(ip,0);

      // No calc needed if swc==0 and no jacobians requested.
      if( (rwc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Negative rwc?
      Numeric psd_weight = 1.0;
      if( rwc < 0 )
        {
          psd_weight = -1.0;
          rwc       *= -1.0;
        }
      
      // Calculate PSD
      Vector psd_1p(nsi);
      if( rwc != 0 )
        {
          psd_rain_W16 ( psd_1p, psd_size_grid, rwc );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the latter value becomes the perturbation for RWC=0.
          const Numeric drwc = max( 0.001*rwc, 1e-7 );
          const Numeric rwcp = rwc + drwc;
          psd_rain_W16 ( psd_1p, psd_size_grid, rwcp );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / drwc; }
        }   
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFromParticleBulkProps(
         Workspace&                   ws,
         Tensor4&                     pnd_field,
         ArrayOfTensor4&              dpnd_field_dx,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Vector&                      lat_grid,
   const Vector&                      lon_grid,
   const Tensor3&                     t_field,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const ArrayOfString&               scat_species,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const ArrayOfArrayOfScatteringMetaData&   scat_meta,
   const Tensor4&                     particle_bulkprop_field,
   const ArrayOfString&               particle_bulkprop_names,
   const ArrayOfAgenda&               pnd_agenda_array,
   const ArrayOfArrayOfString&        pnd_agenda_array_input_names,
   const Index&                       jacobian_do,
   const ArrayOfRetrievalQuantity&    jacobian_quantities,
   const Verbosity&)
{
  // As we allow this method to be called with cloudbox_checkedCalc, it must
  // contain quite a number of checks.
  
  // Number of scattering species
  const Index nss = scat_data.nelem();

  // Checks (not totally complete, but should cover most mistakes)
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, 
                 p_grid, lat_grid, lon_grid );
  chk_atm_field( "particle_bulkprop_field", particle_bulkprop_field,
                 atmosphere_dim, particle_bulkprop_names.nelem(),
                 p_grid, lat_grid, lon_grid );
  // Further checks of *particle_bulkprop_field* below
  if( !cloudbox_on )
    throw runtime_error( "*cloudbox_on* must be true to use this method." );
  if( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error( "Length of *cloudbox_limits* incorrect with respect "
                         "to *atmosphere_dim*." );
  if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
    throw runtime_error( "Invalid data in pressure part of *cloudbox_limits*." );
  if( atmosphere_dim > 1 )
    {
      if( cloudbox_limits[3]<=cloudbox_limits[2] || cloudbox_limits[2]<0 ||
                                           cloudbox_limits[3]>=lat_grid.nelem() )
        throw runtime_error( "Invalid data in latitude part of *cloudbox_limits*." );
      if( atmosphere_dim > 2 )
        {
          if( cloudbox_limits[5]<=cloudbox_limits[4] || cloudbox_limits[4]<0 ||
                                           cloudbox_limits[5]>=lon_grid.nelem() )
            throw runtime_error( "Invalid data in longitude part of *cloudbox_limits*." );
        }
    }
  if( nss < 1 )
    throw runtime_error( "*scat_data* is empty!." );
  if( scat_species.nelem() != nss )
    throw runtime_error( "*scat_data* and *scat_species* are inconsistent in size." );
  if( scat_meta.nelem() != nss )
    throw runtime_error( "*scat_data* and *scat_meta* are inconsistent in size." );
  if( pnd_agenda_array.nelem() != nss )
    throw runtime_error( "*scat_data* and *pnd_agenda_array* are inconsistent "
                         "in size." );
  if( pnd_agenda_array_input_names.nelem() != nss )
    throw runtime_error( "*scat_data* and *pnd_agenda_array_input_names* are "
                         "inconsistent in size." );
  // Further checks of scat_data vs. scat_meta below  

  
  // Effective lengths of cloudbox
  const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index ip_offset = cloudbox_limits[0];
  Index nlat = 1;
  Index ilat_offset = 0;
  if( atmosphere_dim > 1 )
    {
      nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      ilat_offset = cloudbox_limits[2];
    }
  Index nlon = 1;
  Index ilon_offset = 0;
  if( atmosphere_dim > 2 )
    {
      nlat = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      ilon_offset = cloudbox_limits[4];
    }

  // Check that *particle_bulkprop_field* contains zeros outside and at
  // cloudbox boundaries
  const String estring = "*particle_bulkprop_field* can only contain non-zero "
    "values inside the cloudbox.";
  // Pressure end ranges
  for( Index ilon=0; ilon<nlon; ilon++ ) {
    for( Index ilat=0; ilat<nlat; ilat++ ) {
      for( Index ip=0; ip<=cloudbox_limits[0]; ip++ ) {
        if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
          throw runtime_error( estring ); } 
      for( Index ip=cloudbox_limits[1]; ip<p_grid.nelem(); ip++ ) {
        if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
          throw runtime_error( estring ); } } }
  if( atmosphere_dim > 1 )
    {
      // Latitude end ranges
      for( Index ilon=0; ilon<nlon; ilon++ ) {
        for( Index ip=cloudbox_limits[0]+1; ip<cloudbox_limits[1]-1; ip++ ) {
          for( Index ilat=0; ilat<=cloudbox_limits[2]; ilat++ ) {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
              throw runtime_error( estring ); }
          for( Index ilat=cloudbox_limits[3]; ilat<lat_grid.nelem(); ilat++ ) {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
              throw runtime_error( estring ); } } }
      if( atmosphere_dim > 2 )
        {
          // Longitude end ranges
          for( Index ip=cloudbox_limits[0]+1; ip<cloudbox_limits[1]-1; ip++ ) {
            for( Index ilat=cloudbox_limits[2]+1; ilat<cloudbox_limits[3]-1; ilat++ ) {
              for( Index ilon=0; ilon<=cloudbox_limits[4]; ilon++ ) {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                  throw runtime_error( estring ); }
              for( Index ilon=cloudbox_limits[5]; ilon<lon_grid.nelem(); ilon++ ) {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                  throw runtime_error( estring ); } } }
        }
    }
  
  // Cumulative number of scattering elements
  ArrayOfIndex ncumse(nss+1);
  ncumse[0] = 0;
  for( Index i=0; i<nss; i++ )
    {
      if( scat_data[i].nelem() != scat_meta[i].nelem() )
        throw runtime_error( "*scat_data* and *scat_meta* have inconsistent sizes." );
      ncumse[i+1] = ncumse[i] + scat_data[i].nelem();
    }

  // Allocate output variables
  //
  pnd_field.resize( ncumse[nss], np, nlat, nlon );
  pnd_field = 0.0;  // To set all end values to zero
  //
  // Help variables for partial derivatives
  Index                nq = 0;
  ArrayOfArrayOfIndex  scatspecies_to_jq;
  //
  if( !jacobian_do )
    { dpnd_field_dx.resize(0); }
  else
    {
      nq = jacobian_quantities.nelem();
      dpnd_field_dx.resize( nq );
      scatspecies_to_jq.resize( nss );
      //
      for( Index iq=0; iq<nq; iq++ )
        {
          if( jacobian_quantities[iq].MainTag() == SCATSPECIES_MAINTAG  )
            {
              const Index ihit = find_first( scat_species,
                                             jacobian_quantities[iq].Subtag() );
              if( ihit < 0 )
                {
                  ostringstream os;
                  os << "Jacobian quantity with index " << iq << " refers to\n"
                     << "  " << jacobian_quantities[iq].Subtag()
                     << "\nbut this species could not be found in *scat_species*.";
                  throw runtime_error(os.str());
                }
              scatspecies_to_jq[ihit].push_back( iq );
              dpnd_field_dx[iq].resize( ncumse[nss], np, nlat, nlon );
              dpnd_field_dx[iq] = 0.0;  // To set all end values to zero
            }
        }
    }

  // Extract data from pnd-agenda array
  for( Index is=0; is<nss; is++ )
    {
      // Index range with respect to pnd_field
      Range se_range( ncumse[is], ncumse[is+1]-ncumse[is] );

      // Determine how pnd_agenda_array_input_names are related to input fields
      //
      const Index nin = pnd_agenda_array_input_names[is].nelem();
      ArrayOfIndex i_pbulkprop(nin);
      //
      for( Index i=0; i<nin; i++ )
        {
          // We flag temperature with -100
          if( pnd_agenda_array_input_names[is][i] == "Temperature" )
            { i_pbulkprop[i] = -100; }
          else
            {
              i_pbulkprop[i] = find_first( particle_bulkprop_names,
                                           pnd_agenda_array_input_names[is][i] );
              if( i_pbulkprop[i] < 0 )
                {
                  ostringstream os;
                  os << "Pnd-agenda with index " << is << " is set to require \""
                     << pnd_agenda_array_input_names[is][i] << "\",\nbut this quantity "
                     << "could not found in *particle_bulkprop_names*.\n"
                     << "(Note that temperature must be written as \"Temperature\")";
                  throw runtime_error(os.str());
                }
            }
        }

      // Set *dpnd_data_dx_names*
      //
      Index ndx = 0;
      ArrayOfString dpnd_data_dx_names( 0 );
      //
      if( jacobian_do )
        {
          ndx = scatspecies_to_jq[is].nelem();
          dpnd_data_dx_names.resize( ndx );
          for( Index ix=0; ix<ndx; ix++ )
            { dpnd_data_dx_names[ix] =
                jacobian_quantities[scatspecies_to_jq[is][ix]].SubSubtag(); }
        }
      
      // Loop lat/lon positions and call *pnd_agenda*
      for( Index ilon=0; ilon<nlon; ilon++ )
        { 
          for( Index ilat=0; ilat<nlat; ilat++ )
            {
              // Note that we don't need any calculations for end points

              // Here we consider this for lat and lon
              if( ( nlat > 1  &&  ( ilat == 0  ||  ilat == nlat-1 ) )  || 
                  ( nlon > 1  &&  ( ilon == 0  ||  ilon == nlon-1 ) ) )
                { continue; }
                  
              // Pressure handled here, by not including end points in loops
              Matrix pnd_agenda_input( np-2, nin );

              
              for( Index i=0; i<nin; i++ )
                {
                  if( i_pbulkprop[i] == -100 )
                    {
                      for( Index ip=1; ip<np-1; ip++ )
                        { pnd_agenda_input(ip-1 ,i) = t_field(
                                                         ip_offset   + ip,
                                                         ilat_offset + ilat,
                                                         ilon_offset + ilon ); }
                    }
                  else
                    {
                      for( Index ip=1; ip<np-1; ip++ )
                        { pnd_agenda_input(ip-1,i) = particle_bulkprop_field(
                                                         i_pbulkprop[i],
                                                         ip_offset   + ip,
                                                         ilat_offset + ilat,
                                                         ilon_offset + ilon );
                        }
                    }
                }
              
              // Call pnd-agenda array
              Matrix pnd_data;
              Tensor3 dpnd_data_dx;
              //
              pnd_agenda_arrayExecute( ws, pnd_data, dpnd_data_dx, is,
                                       pnd_agenda_input,
                                       pnd_agenda_array_input_names[is],
                                       dpnd_data_dx_names, pnd_agenda_array );

              // Copy to output variables
              for( Index ip=1; ip<np-1; ip++ )
                { pnd_field(se_range,ip,ilat,ilon) = pnd_data(ip-1,joker); }
              for( Index ix=0; ix<ndx; ix++ )
                {
                  for( Index ip=1; ip<np-1; ip++ )
                    { dpnd_field_dx[scatspecies_to_jq[is][ix]]
                            (se_range,ip,ilat,ilon) = dpnd_data_dx(ix,ip-1,joker); }
                }
            }
        }
    }
}
   
                           



/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_W16 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_mass_equivalent,
              const Numeric& RWC,
              const Index& robust,
              const Verbosity&)
{
  Index n_se = diameter_mass_equivalent.nelem();
  dNdD.resize(n_se);

  // abort if RWC is negative
  if ( !robust && RWC<0. )
    {
      ostringstream os;
      os << "Rain water content can not be negative."
         << " Yours is " << RWC << "kg/m3.";
      throw runtime_error ( os.str() );
    }
  Numeric rwc = max( RWC, 0. );

  // calculate particle size distribution with W16
  // [# m^-3 m^-1]
  psd_rain_W16 ( dNdD, diameter_mass_equivalent, rwc );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MH97 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_mass_equivalent,
              const Numeric& IWC,
              const Numeric& T,
              const Index& noisy,
              const Index& robust,
              const Verbosity&)
{
  Index n_se = diameter_mass_equivalent.nelem();
  dNdD.resize(n_se);

  // abort if IWC is negative
  if ( !robust && IWC<0. )
    {
      ostringstream os;
      os << "Ice water content can not be negative."
         << " Yours is " << IWC << "kg/m3.";
      throw runtime_error ( os.str() );
    }
  Numeric iwc = max( IWC, 0. );

  // abort if T is negative
  if ( T<0. )
    {
      ostringstream os;
      os << "Negative temperatures not allowed.\n"
         << "Yours is " << T << "K.";
      throw runtime_error ( os.str() );
    }
  // abort if T is too high
  if ( !robust && T>280. )
    {
      ostringstream os;
      os << "Temperatures above 280K not allowed by MH97"
         << " (to allow: run with robust option).\n"
         << "Yours is " << T << "K.";
      throw runtime_error ( os.str() );
    }
  // allow some margin on T>0C (but use T=0C for PSD calc)
  Numeric t = min( T, 273.15 );
  //Numeric t = T;

  // calculate particle size distribution with MH97
  // [# m^-3 m^-1]
  psd_cloudice_MH97 ( dNdD, diameter_mass_equivalent, iwc, t, noisy );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_F07 (//WS Output:
                Vector& dNdD,
                //WS Input:
                const Vector& diameter_max,
                const Numeric& SWC,
                const Numeric& T,
                const String& regime,
                const Numeric& alpha,
                const Numeric& beta,
                const Index& robust,
                const Verbosity&)
{
  Index n_se = diameter_max.nelem();
  dNdD.resize(n_se);
    
  // abort if SWC is negative
  if ( !robust && SWC<0. )
    {
      ostringstream os;
      os << "Snow water content can not be negative."
         << " Yours is " << SWC << "kg/m3.";
      throw runtime_error ( os.str() );
    }
  Numeric swc = max( SWC, 0. );

  // abort if T is negative
  if ( T<0. )
    {
      ostringstream os;
      os << "Negative temperatures not allowed.\n"
         << "Yours is " << T << "K.";
      throw runtime_error ( os.str() );
    }

  // check regime
  if( regime!="TR" && regime!="ML" )
    {
      ostringstream os;
      os << "Only regimes \"TR\" (tropical) and \"ML\" (midlatitude) known.\n"
         << "Yours is \"" << regime << "\".";
      throw runtime_error ( os.str() );
    }

  // calculate particle size distribution with F07TR
  // [# m^-3 m^-1]
  psd_snow_F07(dNdD, diameter_max, swc, T, alpha, beta, regime);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H11 (//WS Output:
             Vector& dNdD,
             //WS Input:
             const Vector& Dmax,
             const Numeric& t,
             const Verbosity&)
{
  Index n_se = Dmax.nelem();
  dNdD.resize(n_se);

  for ( Index i=0; i<n_se; i++ ) //loop over number of scattering elementss
    {
      // calculate particle size distribution for H11
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_H11 ( Dmax[i], t );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_Ar_H13 (//WS Output:
                Vector& dNdD,
                Vector& Ar,
                //WS Input:
                const Vector& Dmax,
                const Numeric& t,
                const Verbosity&)
{
  Index n_se = Dmax.nelem();
  dNdD.resize(n_se);
  Ar.resize(n_se);

  for ( Index i=0; i<n_se; i++ ) //loop over number of scattering elementss
    {
      // calculate particle size distribution for H13Shape
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_H13Shape ( Dmax[i], t );
      // calculate Area ratio distribution for H13Shape
      Ar[i] = area_ratioH13 ( Dmax[i], t );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_S2M (//WS Output:
                 Vector& dNdD,
                 //WS Input:
               const Vector& mass,
               const Numeric& N_tot,
               const Numeric& M,
               const String& psd_type,
               const Verbosity&)
{
    Index n_se = mass.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with S2M
        // [# m^-3 kg^-1]
        dNdD[i] = WCtopnd_S2M(mass[i], N_tot, M, psd_type) ;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_S2M_M (//WS Output:
               Vector& dNdD,
               //WS Input:
               const Vector& mass,
               const Numeric& mean_mass,
               const Numeric& M,
               const String& psd_type,
               const Verbosity&)
{
    Numeric N_tot;

    // Calculate the total number density from mass density M and the
    // mean particle mass
    N_tot=M/mean_mass;
    
    Index n_se = mass.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with S2M
        // [# m^-3 kg^-1]
        dNdD[i] = WCtopnd_S2M(mass[i], N_tot, M, psd_type) ;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MGD_LWC (//WS Output:
                 Vector& dNdD,
                 //WS Input:
                 const Vector& diameter_volume_equ,
                 const Numeric& rho,
                 const Numeric& LWC,
                 const Verbosity&)
{
    Index n_se = diameter_volume_equ.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with ModGamma for liquid
        // [# m^-3 m^-1]
        dNdD[i] = LWCtopnd_MGD_LWC( diameter_volume_equ[i],rho ,LWC );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MGD_IWC (//WS Output:
                   Vector& dNdD,
                   //WS Input:
                   const Vector& diameter_volume_equ,
                   const Numeric& rho,
                   const Numeric& IWC,
                   const Verbosity&)
{
    Index n_se = diameter_volume_equ.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with ModGamma for ice
        // [# m^-3 m^-1]
        dNdD[i] = IWCtopnd_MGD_IWC( diameter_volume_equ[i],rho,IWC );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H98 (//WS Output:
             Vector& dNdD,
             //WS Input:
             const Vector& diameter_volume_equivalent,
             const Numeric& LWC,
             const Verbosity&)
{
  Index n_se = diameter_volume_equivalent.nelem();
  dNdD.resize(n_se);
  const Numeric dDdR = 2.; 

  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution for liquid
      // and compensate for LWCtopnd providing dNdR
      // [# m^-3 m^-1]
      dNdD[i] = LWCtopnd ( LWC, diameter_volume_equivalent[i]/2. ) / dDdR;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MP48 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_melted_equivalent,
              const Numeric& PR,
              const String& PRunit,
              const Numeric& rho,
              const Verbosity&)
{
  Numeric tPR;
  if (PRunit == "mm/h")
    tPR = PR;
  else if ( (PRunit == "SI") || (PRunit == "kg/m2/s") )
    {
      if (rho<=0.)
        {
          ostringstream os;
          os << "Precipitation unit " << PRunit
             << " requires valid material density (rho>0).\n"
             << "Yours is rho=" << rho << "kg/m3.\n";
          throw runtime_error ( os.str() );
        }
      tPR = PR * (3.6e6/rho);
    }
  else
    {
      ostringstream os;
      os << "Precipitation unit '" << PRunit << "' unknown.\n";
      throw runtime_error ( os.str() );
    }

  Index n_se = diameter_melted_equivalent.nelem();
  dNdD.resize(n_se);

  // derive particle number density for all given sizes
  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution with MP48
      // output: [# m^-3 m^-1]
      dNdD[i] = PRtopnd_MP48 ( tPR, diameter_melted_equivalent[i]);
    }
}


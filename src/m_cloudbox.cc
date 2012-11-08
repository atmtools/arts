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
#include "sorting.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Numeric DEG2RAD;
extern const Numeric PI;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxOff (
         Index&                       cloudbox_on,
         ArrayOfIndex&                cloudbox_limits,
         Agenda&                      iy_cloudbox_agenda,
         Tensor4&                     pnd_field,
         Index&                       use_mean_scat_data,
         ArrayOfSingleScatteringData& scat_data_raw,
         Matrix&                      particle_masses,
   const Verbosity&)
{
  cloudbox_on = 0;
  cloudbox_limits.resize ( 0 );
  iy_cloudbox_agenda = Agenda();
  iy_cloudbox_agenda.set_name ( "iy_cloudbox_agenda" );
  pnd_field.resize(0,0,0,0);
  use_mean_scat_data = -999;
  scat_data_raw.resize(0);
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
                               const Tensor4&  massdensity_field,
                               // Control Parameters
                               const Numeric&  cloudbox_margin,
                               const Verbosity& verbosity)
{
  // Variables
  //const Index np = massdensity_field.npages(); 
  Index p1 = massdensity_field.npages()-1;
  Index p2 = 0;
  DEBUG_ONLY(
             Index lat1 = massdensity_field.nrows()-1;
             Index lat2 = 0;
             Index lon1 = massdensity_field.ncols()-1;
             Index lon2 = 0;
  )

  Index i, j, k, l;

  // initialize flag, telling if all selected *massdensity_fields* are
  // zero(false) or not(true)
  bool x = false; 


  // Check existing WSV
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  // includes p_grid chk_if_decresasing
  chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid ); 
  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize ( atmosphere_dim*2 );

  //--------- Start loop over particles ----------------------------------------
  for ( l=0; l<massdensity_field.nbooks(); l++ )
  {
    bool y; // empty_flag
    //y is set to true, if a single value of massdensity_field is unequal zero
    chk_massdensity_field ( y,
                         atmosphere_dim,
                         massdensity_field ( l, joker, joker, joker ),
                         p_grid,
                         lat_grid,
                         lon_grid );

    //-----------Start setting cloudbox limits----------------------------------
    if ( y )
    {
      //massdensity_field unequal zero -> x is true
      x = true;

      if ( atmosphere_dim == 1 )
      {
        // Pressure limits
        ConstVectorView hydro_p = massdensity_field ( l, joker, 0 , 0 );


        // set lower cloudbox_limit to surface if margin = -1
        if ( cloudbox_margin == -1 )
        {
          cloudbox_limits[0] = 0;
          i = p1 = 0;
        }
        else
        {
          // find index of first pressure level where hydromet_field is
          // unequal 0, starting from the surface
          for ( i=0; i<hydro_p.nelem(); i++ )
          {
            if ( hydro_p[i] != 0.0 )
            {
               // check if p1 is the lowest index in all selected
               // massdensity fields
               if ( p1 > i )
              {
                p1 = i;
              }
              break;
            }
          }

        }
        // find index of highest pressure level, where massdensity_field is
        // unequal 0, starting from top of the atmosphere
        for ( j=hydro_p.nelem()-1; j>=i; j-- )
        {
          if ( hydro_p[j] != 0.0 )
          {
             // check if p2 is the highest index in all selected
             // massdensity fields
             if ( p2 < j )
            {
	      p2 = j;
            }
            break;
          }
        }

      }
    }

    /*  //NOT WORKING YET
      // Latitude limits
      else if ( atmosphere_dim == 2 )
      {
        MatrixView hydro_lat = hydromet_field ( nhyd, joker, joker, 0 );

        for ( i=0; i<hydro_lat.nrows(); i++ )
        {
          for ( j=0; j<hydro_lat.ncols(); j++ )
          {
            if ( hydro_lat[i,j] != 0.0 )
            {

              if ( lat1 <= j ) lat1 =j;
              //cloudbox_limits[2] = lat1;
              //break;
            }

          }
          if ( p1 <= i )    p1 = i;
        }

        for ( k=hydro_lat.nelem()-1; k>=i; k-- )
        {
          if ( hydro_lat[k] != 0.0 )
          {
            lat2 = k;
            cloudbox_limits[3] = lat2;
            break;

          }

        }
      }

      // Longitude limits
      if ( atmosphere_dim == 3 )
      {
        Tensor3View hydro_lon = hydromet_field ( nhyd, joker, joker, joker );

        for ( i=0; i<hydro_lon.nelem(); i++ )
        {
          if ( hydro_lon[i] != 0.0 )
          {
            lon1 = i;
            cloudbox_limits[4] = lon1;
            break;
          }

        }
        for ( j=hydro_lon.nelem()-1; j>=i; j-- )
        {
          if ( hydro_lon[j] != 0.0 )
          {
            lon2 = j;
            cloudbox_limits[5] = lon2;
            break;

          }


}*/
  }
  // decrease lower cb limit by one to ensure that linear interpolation of 
  // particle number densities is possible.
  Index p0 = 0; //only for the use of function *max*
  p1 = max(p1-1, p0);

  Numeric p_margin1;

  // alter lower cloudbox_limit by cloudbox_margin, using barometric
  // height formula
  p_margin1 = barometric_heightformula ( p_grid[p1], cloudbox_margin );
  k = 0;
  while ( p_grid[k+1] >= p_margin1 && k+1 < p_grid.nelem() ) k++;
  cloudbox_limits[0]= k;

  // increase upper cb limit by one to ensure that linear interpolation of 
  // particle number densities is possible.
  p2 = min(p2+1, massdensity_field.npages()-1);
  // set upper cloudbox_limit
  // if cloudbox reaches to the upper most pressure level
  if ( p2 >= massdensity_field.npages()-1)
  {
    CREATE_OUT2;
    out2<<"The cloud reaches to TOA!\n"
    <<"Check massdensity_field data, if realistic!\n";
  }
  cloudbox_limits[1] = p2;

  //out0<<"\n"<<p2<<"\n"<<p_grid[p2]<<"\n";

  // check if all selected massdensity fields are zero at each level,
  // than switch cloudbox off, skipping scattering calculations
  if ( !x )
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
void cloudbox_checkedCalc(      
         Index&          cloudbox_checked,
   const Index&          basics_checked,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        z_field,
   const Matrix&         z_surface,
   const Tensor3&        wind_u_field,
   const Tensor3&        wind_v_field,
   const Tensor3&        wind_w_field,
   const Index&          cloudbox_on,    
   const ArrayOfIndex&   cloudbox_limits,
   const Tensor4&        pnd_field,
   const ArrayOfSingleScatteringData& scat_data_raw,
   const Matrix&         particle_masses,
   const Verbosity&)
{
  // Demanded space between cloudbox and lat and lon edges [degrees]
  const Numeric llmin = 20;

  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  
  chk_if_bool( "cloudbox_on", cloudbox_on );

  if( cloudbox_on )
    {
      // Winds, must be empty variables (i.e. no winds allowed)
      {
        ostringstream ow;
        ow << "The scattering methods are not (yet?) handling winds. For this\n"
           << "reason, the WSVs for wind fields must all be empty with an\n."
           << "active cloudbox.";
        if( wind_w_field.npages() > 0 )
          { throw runtime_error( ow.str() ); }
        if( wind_v_field.npages() > 0 )
          { throw runtime_error( ow.str() ); }
        if( atmosphere_dim > 2  &&  wind_u_field.npages() > 0 )
          { throw runtime_error( ow.str() ); }
      }

      // Cloudbox limits
      if( cloudbox_limits.nelem() != atmosphere_dim*2 )
        {
          ostringstream os;
          os << "The array *cloudbox_limits* has incorrect length.\n"
             << "For atmospheric dim. = " << atmosphere_dim 
             << " the length shall be " << atmosphere_dim*2
             << " but it is " << cloudbox_limits.nelem() << ".";
          throw runtime_error( os.str() );
        }
       if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
        {
          ostringstream os;
          os << "Incorrect value(s) for cloud box pressure limit(s) found."
             << "\nValues are either out of range or upper limit is not "
             << "greater than lower limit.\nWith present length of "
             << "*p_grid*, OK values are 0 - " << p_grid.nelem()-1
             << ".\nThe pressure index limits are set to " 
             << cloudbox_limits[0] << " - " << cloudbox_limits[1] << ".";
          throw runtime_error( os.str() );
        }
      if( atmosphere_dim >= 2 )
        {
          const Index n = lat_grid.nelem();
          if( cloudbox_limits[3]<=cloudbox_limits[2] || cloudbox_limits[2]<1 ||
                                                      cloudbox_limits[3]>=n-1 )
            {
              ostringstream os;
              os << "Incorrect value(s) for cloud box latitude limit(s) found."
                 << "\nValues are either out of range or upper limit is not "
                 << "greater than lower limit.\nWith present length of "
                 << "*lat_grid*, OK values are 1 - " << n-2
                 << ".\nThe latitude index limits are set to " 
                 << cloudbox_limits[2] << " - " << cloudbox_limits[3] << ".";
              throw runtime_error( os.str() );
            }
          if( ( lat_grid[cloudbox_limits[2]] - lat_grid[0] < llmin )  &&
              ( atmosphere_dim==2  ||  
              ( atmosphere_dim==3 && lat_grid[0]>-90) ) )
            {
              ostringstream os;
              os << "Too small distance between cloudbox and lower end of "
                 << "latitude grid.\n"
                 << "This distance must be " << llmin << " degrees.\n"
                 << "Cloudbox ends at " << lat_grid[cloudbox_limits[2]]
                 << " and latitude grid starts at " << lat_grid[0] << ".";
              throw runtime_error( os.str() );
            }
          if( ( lat_grid[n-1] - lat_grid[cloudbox_limits[3]] < llmin )  &&
              ( atmosphere_dim==2  ||  
              (atmosphere_dim==3 && lat_grid[n-1]<90) ) )
            {
              ostringstream os;
              os << "Too small distance between cloudbox and upper end of "
                 << "latitude grid.\n"
                 << "This distance must be " << llmin << " degrees.\n"
                 << "Cloudbox ends at " << lat_grid[cloudbox_limits[3]]
                 << " and latitude grid ends at " << lat_grid[n-1] << ".";
              throw runtime_error( os.str() );
            }
        }
      if( atmosphere_dim >= 3 )
        {
          const Index n = lon_grid.nelem();
          if( cloudbox_limits[5]<=cloudbox_limits[4] || cloudbox_limits[4]<1 ||
                                                      cloudbox_limits[5]>=n-1 )
            {
              ostringstream os;
              os << "Incorrect value(s) for cloud box longitude limit(s) found"
                 << ".\nValues are either out of range or upper limit is not "
                 << "greater than lower limit.\nWith present length of "
                 << "*lon_grid*, OK values are 1 - " << n-2
                 << ".\nThe longitude limits are set to " 
                 << cloudbox_limits[4] << " - " << cloudbox_limits[5] << ".";
              throw runtime_error( os.str() );
            }
          if( lon_grid[n-1] - lon_grid[0] < 360 )
            {
              const Numeric latmax = max( abs(lat_grid[cloudbox_limits[2]]),
                                          abs(lat_grid[cloudbox_limits[3]]) );
              const Numeric lfac = 1 / cos( DEG2RAD*latmax );
              if( lon_grid[cloudbox_limits[4]]-lon_grid[0] < llmin/lfac )
                {
                  ostringstream os;
                  os << "Too small distance between cloudbox and lower end of\n"
                     << "longitude grid. This distance must here be " 
                     << llmin/lfac << " degrees.";
                  throw runtime_error( os.str() );
                }
              if( lon_grid[n-1] - lon_grid[cloudbox_limits[5]] < llmin/lfac )
                {
                  ostringstream os;
                  os << "Too small distance between cloudbox and upper end of\n"
                     << "longitude grid. This distance must here be " 
                     << llmin/lfac << " degrees.";
                  throw runtime_error( os.str() );
                }
            }
        }

      // pnd_field
      //
      const Index np = scat_data_raw.nelem();
      // Dummy variables to mimic grids of correct size
      Vector g1( cloudbox_limits[1]-cloudbox_limits[0]+1 ), g2(0), g3(0);
      if( atmosphere_dim >= 2 ) 
        { g2.resize( cloudbox_limits[3]-cloudbox_limits[2]+1 ); }
      if( atmosphere_dim == 3 ) 
        { g3.resize( cloudbox_limits[5]-cloudbox_limits[4]+1 ); }
      //
      chk_atm_field( "pnd_field", pnd_field, atmosphere_dim, np, g1, g2, g3 );
      //
      if( min(pnd_field) < 0 )
        throw runtime_error( "All values in *pnd_field* must be >= 0." );
      //
      for( Index a=0; a<g2.nelem(); a++ ) { 
        for( Index o=0; o<g3.nelem(); o++ ) { 
          if( max(pnd_field(joker,0,a,o)) > 0  && 
              z_field(cloudbox_limits[0],a,o) > z_surface(a,o) )
            throw runtime_error( "A non-zero value found in *pnd_field* at the "
                             "lower altitude limit of the cloudbox (but the "
                             "position is not below the surface altitude)." );
        } }
      if( max(pnd_field(joker,g1.nelem()-1,joker,joker)) > 0 )
        throw runtime_error( "A non-zero value found in *pnd_field* at "
                             "upper altitude limit of the cloudbox." );
      if( atmosphere_dim >= 2 )
        {
          if( max(pnd_field(joker,joker,0,joker)) > 0 )
            throw runtime_error( "A non-zero value found in *pnd_field* at "
                                 "lower latitude limit of the cloudbox." );
          if( max(pnd_field(joker,joker,g2.nelem()-1,joker)) > 0 ) 
            throw runtime_error( "A non-zero value found in *pnd_field* at "
                                 "upper latitude limit of the cloudbox." );
        }
      if( atmosphere_dim == 3 )
        {
          if( max(pnd_field(joker,joker,joker,0)) > 0 )
            throw runtime_error( "A non-zero value found in *pnd_field* at "
                                 "lower longitude limit of the cloudbox." );
          if( max(pnd_field(joker,joker,joker,g3.nelem()-1)) > 0 ) 
            throw runtime_error( "A non-zero value found in *pnd_field* at "
                                 "upper longitude limit of the cloudbox." );
        }

      // particle_masses
      //
      if( particle_masses.nrows() > 0 )
        {
          if( particle_masses.nrows() != np )
            throw runtime_error( "The WSV *particle_masses* must either be "
                                 "empty or have a row size matching the "
                                 "length of *scat_data_raw*." );
          if( min(particle_masses) < 0 )
            throw runtime_error( 
                            "All values in *particles_masses* must be >= 0." );
        }
    }

  // If here, all OK
  cloudbox_checked = 1;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void Massdensity_cleanup (//WS Output:
                          Tensor4& massdensity_field,
                          //WS Input:
                          const Numeric& massdensity_threshold,
                          const Verbosity&)
{
  // Check that massdensity_fields contain realistic values
  //(values smaller than massdensity_threshold will be set to 0)
  for ( Index i=0; i<massdensity_field.nbooks(); i++ )
  {
    for ( Index j=0; j<massdensity_field.npages(); j++ )
    {
      for ( Index k=0; k<massdensity_field.nrows(); k++ )
      {
        for ( Index l=0; l<massdensity_field.ncols(); l++ )
        {
          if ( massdensity_field ( i,j,k,l ) < massdensity_threshold ) 
          {
            massdensity_field ( i,j,k,l ) = 0.0;
          }
        }
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleSpeciesInit (ArrayOfString&  part_species,
                          const Verbosity&)
{
  part_species.resize ( 0 );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleSpeciesSet (// WS Generic Output:
                         ArrayOfString&  part_species,
                         // Control Parameters:
                         const ArrayOfString& particle_tags,
                         const String& delim,
                         const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  part_species.resize ( particle_tags.nelem() );
  //assign input strings to part_species
  part_species = particle_tags;

  chk_part_species (part_species, delim);

  // Print list of particle settings to the most verbose output stream:
  out3 << "  Defined particle settings: ";
  for ( Index i=0; i<part_species.nelem(); ++i )
  {
    out3 << "\n  " << i << ": "<<part_species[i];

  }
  out3 << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeInit (//WS Output:
                       ArrayOfSingleScatteringData& scat_data_raw,
                       ArrayOfGriddedField3& pnd_field_raw,
                       const Verbosity&)
{
  scat_data_raw.reserve ( 20 );
  pnd_field_raw.reserve ( 20 );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAddAll (//WS Output:
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
                         const String& pnd_field_file,
                         const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Cloudbox limits
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error (
      "*cloudbox_limits* is a vector which contains"
      "the upper and lower\n"
      "limit of the cloud for all atmospheric dimensions.\n"
      "So its length must be 2 x *atmosphere_dim*" );
  // Frequency grid
  if ( f_grid.nelem() == 0 )
    throw runtime_error ( "The frequency grid is empty." );
  chk_if_increasing ( "f_grid", f_grid );


  //--- Reading the data ---------------------------------------------------
  ArrayOfString data_files;
  xml_read_from_file ( filename_scat_data, data_files, verbosity );
  scat_data_raw.resize ( data_files.nelem() );

  for ( Index i = 0; i<data_files.nelem(); i++ )
  {

    out2 << "  Read single scattering data\n";
    xml_read_from_file ( data_files[i], scat_data_raw[i], verbosity );

    chk_single_scattering_data ( scat_data_raw[i],
                                 data_files[i], f_grid,
                                 verbosity );

  }

  out2 << "  Read particle number density data \n";
  xml_read_from_file ( pnd_field_file, pnd_field_raw, verbosity );

  chk_pnd_raw_data ( pnd_field_raw,
                     pnd_field_file, atmosphere_dim, p_grid, lat_grid,
                     lon_grid, cloudbox_limits, verbosity);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringParticleTypeAndMetaRead (//WS Output:
                                        ArrayOfSingleScatteringData& scat_data_raw,
                                        ArrayOfScatteringMetaData& scat_data_meta_array,
                                        const Vector& f_grid,
                                        // Keywords:
                                        const String& filename_scat_data,
                                        const String& filename_scat_meta_data,
                                        const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  //--- Reading the data ---------------------------------------------------
  ArrayOfString data_files;
  ArrayOfString meta_data_files;
  
  // single scattering data read to temporary ArrayOfSingleScatteringData
  xml_read_from_file ( filename_scat_data, data_files, verbosity );
  scat_data_raw.resize ( data_files.nelem() );

  for ( Index i = 0; i<data_files.nelem(); i++ )
  {
    out3 << "  Read single scattering data\n";
    xml_read_from_file ( data_files[i], scat_data_raw[i], verbosity );

    chk_single_scattering_data ( scat_data_raw[i],
                                 data_files[i], f_grid,
                                 verbosity );

  }

  // scattering meta data read to temporary ArrayOfScatteringMetaData
  xml_read_from_file ( filename_scat_meta_data, meta_data_files, verbosity );
  scat_data_meta_array.resize ( meta_data_files.nelem() );

  for ( Index i = 0; i<meta_data_files.nelem(); i++ )
  {

    out3 << "  Read scattering meta data\n";
    xml_read_from_file ( meta_data_files[i], scat_data_meta_array[i], verbosity );

    //FIXME: currently nothing is done in chk_scattering_meta_data!
    chk_scattering_meta_data ( scat_data_meta_array[i],
                               meta_data_files[i], verbosity );

  }
  
  // check if arrays have same size
  chk_scattering_data ( scat_data_raw,
                        scat_data_meta_array, verbosity );
  
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringParticlesSelect (//WS Output:
                                ArrayOfSingleScatteringData& scat_data_raw,
                                ArrayOfScatteringMetaData& scat_data_meta_array,
                                ArrayOfIndex& scat_data_nelem,
                                // WS Input:
                                const ArrayOfString& part_species,
                                const String& delim,
                                const Verbosity& verbosity)
{ 
  CREATE_OUT1;
  CREATE_OUT3;
  //--- Adjusting data to user specified input (part_species)-------------------
  
  Index intarr_total = 0;
  ArrayOfIndex intarr;
  
  // make temporary copy
  ArrayOfSingleScatteringData scat_data_raw_tmp = scat_data_raw;
  ArrayOfScatteringMetaData scat_data_meta_array_tmp = scat_data_meta_array;
  
  scat_data_nelem.resize( part_species.nelem() );
  
  ArrayOfIndex selected;
  selected.resize(scat_data_meta_array_tmp.nelem());
  selected = 0;
  // loop over array of part_species--------------------------------------------
  for ( Index k=0; k<part_species.nelem(); k++ )
  {
   
    String part_type;
    String partfield_name;
    Numeric sizemin;
    Numeric sizemax;

    //split part_species string and copy values to parameter
    parse_partfield_name( partfield_name, part_species[k], delim);
    parse_part_type( part_type, part_species[k], delim);
    parse_part_size(sizemin, sizemax, part_species[k], delim);

    // choosing the specified SingleScatteringData and ScatteringMetaData
    for ( Index j=0; j<scat_data_meta_array_tmp.nelem(); j++ )
    {
      // check for particle (material/phase) type (e.g. "Ice", "Water",...)
      if ( scat_data_meta_array_tmp[j].type == part_type ) 
      {       
        // particle radius is calculated from particle volume given in
        // scattering meta data
        Numeric r_particle = 
          pow ( 3./4. * scat_data_meta_array_tmp[j].V * 1e18 / PI , 1./3. );
	
	// check if particle is in size range
  // (sizemax < 0 results from wildcard usage and means consider all sizes on
  // the upper end)
        if ( r_particle  >= sizemin &&
             ( sizemax >= r_particle || sizemax < 0. ))
	{
	  // fill ArrayOfIndex with indices of selected scattering data
          intarr.push_back ( j );
        }
      selected[j] = 1;
      out3 << "Selecting particle " << j+1 << "/" << scat_data_meta_array_tmp.nelem()
           << " (" << scat_data_meta_array_tmp[j].type << ")\n";
      }
    }
    // WSV scat_data_nelem gets the number of elements of scattering data
    // connected to each selection String in *part_species*   
    scat_data_nelem[k] = intarr.nelem() - intarr_total;
    intarr_total = intarr.nelem();

    // Use particle profile without corresponding ScattData poses high risk of
    // accidentially neglecting those particles. That's unlikely what the user
    // intends. Hence throw error.
    if (scat_data_nelem[k]<1)
      {
        ostringstream os;
        os << "Particle species " << partfield_name << " of type " << part_type
           << " requested.\n"
           << "But no Scattering Data found for it!\n";
        throw runtime_error ( os.str() );
      }
  }

  // check if array is empty
  if ( !intarr.nelem() )
  {
    /* don't throw error. just a warning (should only apply when part_species is empty)!
    ostringstream os;
    os << "The selection in " << part_species << 
        " is NOT choosing any of the given Scattering Data.\n"
        "--> Does the selection in *part_species* fit any of the "
        "Single Scattering Data input? \n";
    throw runtime_error ( os.str() );*/
    out1 << "WARNING! No ScatteringData selected.\n"
         << "Continuing with no particles.\n";
  }

  // check if we ignored any ScatteringMetaData entry
  /* JM120401: what's the meaning of this? does it work properly?
               uncomment again, when answered.
  for ( Index j = 0; j<selected.nelem(); j++)
  {
    if (selected[j]==0)
    {
      out1 << "WARNING! Ignored ScatteringMetaData[" << j << "] ("
           << scat_data_meta_array_tmp[j].type << ")!\n";
    }
  } */


  // resize WSVs to size of intarr
  scat_data_raw.resize ( intarr.nelem() );
  scat_data_meta_array.resize ( intarr.nelem() );

  for ( Index j=0; j<intarr.nelem(); j++ )
  {
    //append to WSV Arrays
    scat_data_meta_array[j] = scat_data_meta_array_tmp[intarr[j]] ;
    scat_data_raw[j] = scat_data_raw_tmp[intarr[j]] ;
  }


}


/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesSet (//WS Output:
                         Matrix& particle_masses,
                         // WS Input:
                         const ArrayOfScatteringMetaData& scat_data_meta_array,
                         const ArrayOfIndex& scat_data_nelem,
                         const ArrayOfString& part_species,
                         const Verbosity& verbosity)
{
  // checks
  if (scat_data_nelem.nelem() != part_species.nelem())
  {
    ostringstream os;
    os << "Dimensions of part_species and scat_data_nelem do not agree.";
    throw runtime_error(os.str());
  }

  // resize particle_masses to required diemsions and properly initialize values
  particle_masses.resize ( part_species.nelem(), scat_data_meta_array.nelem() );
  particle_masses = 0.;
  Index scat_data_start = 0;

  // calculate and set particle_masses
  for ( Index k=0; k<part_species.nelem(); k++ )
  {
    for ( Index j=scat_data_start; j<scat_data_start+scat_data_nelem[k]; j++ )
    {
      particle_masses (k,j) =
        scat_data_meta_array[j].density * scat_data_meta_array[j].V;
    }

    scat_data_start += scat_data_nelem[k];

  }
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
                     const String& pnd_field_file,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  
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
  xml_read_from_file(scat_data_file, scat_data_raw[scat_data_raw.nelem()-1],
                     verbosity);

  chk_single_scattering_data(scat_data_raw[scat_data_raw.nelem()-1],
                             scat_data_file, f_grid, verbosity);       
  
  out2 << "  Read particle number density field\n";
  if (pnd_field_file.nelem() < 1)
  {
    CREATE_OUT1;
    out1 << "Warning: No pnd_field_file specified. Ignored. \n";
  }
  else
    {
      xml_read_from_file(pnd_field_file, pnd_field_raw[pnd_field_raw.nelem()-1],
                         verbosity);
      
      chk_pnd_data(pnd_field_raw[pnd_field_raw.nelem()-1],
                   pnd_field_file, atmosphere_dim, p_grid, lat_grid,
                   lon_grid, cloudbox_limits, verbosity);
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
                   const ArrayOfIndex& cloudbox_limits,
                   const Verbosity&)
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
  const Index Np_cloud = cloudbox_limits[1]-cloudbox_limits[0]+1;
  
  ConstVectorView p_grid_cloud = p_grid[Range(cloudbox_limits[0], Np_cloud)];

  // Check that no scatterers exist outside the cloudbox
  chk_pnd_field_raw_only_in_cloudbox(atmosphere_dim, pnd_field_raw,
                                     p_grid, lat_grid, lon_grid,
                                     cloudbox_limits);
      
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
                  pnd_field_raw[i].data(joker,0,0), gp_p );
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
                  pnd_field_raw[i].data(joker,joker,0), gp_p, gp_lat );
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
                  pnd_field_raw[i].data, gp_p, gp_lat, gp_lon );
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldExpand1D(Tensor4&        pnd_field,
                       const Index&    atmosphere_dim,
                       const Index&    cloudbox_checked,    
                       const Index&    cloudbox_on,    
                       const ArrayOfIndex&   cloudbox_limits,
                       const Index&    nzero,
                       const Verbosity&)
{
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  if( atmosphere_dim == 1 )
    { throw runtime_error( "No use in calling this method for 1D." ); }
  if( !cloudbox_on )
    { throw runtime_error( 
                "No use in calling this method with cloudbox off." ); }

  if( nzero < 1 )
    { throw runtime_error( "The argument *nzero must be > 0." ); }

  // Sizes
  const Index   npart = pnd_field.nbooks();
  const Index   np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index   nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
        Index   nlon = 1;
  if( atmosphere_dim == 3 )  
    { nlon = cloudbox_limits[5] - cloudbox_limits[5] + 1; }

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
                   ArrayOfSingleScatteringData& scat_data_raw,
                   //WS Input:
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const Verbosity&)
{
  // 3D  atmosphere
  if (lat_grid.nelem()>1)
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
  scat_data_raw[0].ptype = PARTICLE_TYPE_MACROS_ISO;
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
void pnd_fieldSetup (//WS Output:
                     Tensor4& pnd_field,
                     //WS Input:
                     const Index& atmosphere_dim,
                     const Index& cloudbox_on,
                     const ArrayOfIndex& cloudbox_limits,
                     const Tensor4& massdensity_field,
                     const Tensor3& t_field,
                     const ArrayOfScatteringMetaData& scat_data_meta_array,
                     const ArrayOfString& part_species,
                     const ArrayOfIndex& scat_data_nelem,
                     const String& delim,
                     const Verbosity& verbosity)
{
  CREATE_OUT1;

  // Cloudbox on/off?
  if ( !cloudbox_on )
  {
    /* Must initialise pnd_field anyway; but empty */
    pnd_field.resize(0, 0, 0, 0);
    return;
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

  if ((limits[1] > massdensity_field.npages()) ||
      (limits[1] > t_field.npages()) ||
      (limits[3] > massdensity_field.nrows()) ||
      (limits[3] > t_field.nrows()) ||
      (limits[5] > massdensity_field.ncols()) ||
      (limits[5] > t_field.ncols()))
  {
    ostringstream os;
    os << "Cloudbox out of bounds compared to fields. "
       << "Upper limits: (p, lat, lon): "
       << "(" << limits[1] << ", " << limits[3] << ", " << limits[5] << "). "
       << "*massdensity_field*: "
       << "(" << massdensity_field.npages() << ", "
       << massdensity_field.nrows() << ", "
       << massdensity_field.ncols() << "). "
       << "*t_field*: "
       << "(" << t_field.npages() << ", "
       << t_field.nrows() << ", "
       << t_field.ncols() << ").";
    throw runtime_error(os.str());
  }

  //resize pnd_field to required atmospheric dimension and scatt particles
  pnd_field.resize ( scat_data_meta_array.nelem(), limits[1]-limits[0],
                     limits[3]-limits[2], limits[5]-limits[4] );
  Index scat_data_start = 0;
  ArrayOfIndex intarr;

  //-------- Start pnd_field calculations---------------------------------------

  // loop over nelem of part_species
  for ( Index k=0; k<part_species.nelem(); k++ )
  {

    String psd_param;
    String partfield_name;
    String part_type;
    String psd;

    //split String and copy to ArrayOfString
    parse_psd_param( psd_param, part_species[k], delim);
    parse_partfield_name( partfield_name, part_species[k], delim);
    parse_part_type( part_type, part_species[k], delim);

    // initialize control parameters
    Vector vol_unsorted ( scat_data_nelem[k], 0.0 );
    Vector d_max_unsorted (scat_data_nelem[k], 0.0);
    Vector vol ( scat_data_nelem[k], 0.0 );
    Vector dm ( scat_data_nelem[k], 0.0 );
    Vector r ( scat_data_nelem[k], 0.0 );
    Vector rho ( scat_data_nelem[k], 0.0 );
    Vector pnd ( scat_data_nelem[k], 0.0 );
    //Vector pnd2 ( scat_data_nelem[k], 0.0 ); //temporary
    Vector dN ( scat_data_nelem[k], 0.0 );
    //Vector dN2 ( scat_data_nelem[k], 0.0 ); //temporary
    //Vector dlwc ( scat_data_nelem[k], 0.0 ); //temporary

    //Tensor3 md_field =  massdensity_field ( k, joker, joker, joker );

    //---- pnd_field calculations for MH97 -------------------------------
    if ( psd_param.substr(0,4) == "MH97" )
    {
        psd = "MH97";

        //check for expected particle field name
        if ( partfield_name != "IWC" )
        {
            out1 << "WARNING! The particle field name is unequal 'IWC'.\n"
                 << psd << " should should only be applied to cloud"
                 << " ice.\n";
        }

        //check for expected particle phase
        if ( part_type != "Ice" )
        {
            out1 << "WARNING! The particle phase is unequal 'Ice'.\n"
                 << psd << " should only be applied to ice"
                 << " particles.\n";
        }

        pnd_fieldMH97( pnd_field,
                       massdensity_field ( k, joker, joker, joker ),
                       t_field, limits,
                       scat_data_meta_array, scat_data_start,
                       scat_data_nelem[k], part_species[k], delim,
                       verbosity);
    }

    //---- pnd_field calculations for H11 ----------------------------
    else if ( psd_param == "H11" )
    {
        psd = "H11";

        //check for expected particle field name
        if ( partfield_name != "IWC" && partfield_name != "Snow")
        {
            out1 << "WARNING! The particle field name is unequal 'IWC' and 'Snow'.\n"
                 << psd << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        //check for expected particle phase
      //check for correct particle phase
        if ( part_type != "Ice" &&  part_type != "Water" )
        {
            out1 << "WARNING! The particle phase is unequal 'Ice'.\n"
                 << psd << " should only be applied to ice (cloud ice or snow)"
                 << " particles.\n";
        }

        pnd_fieldH11 (pnd_field,
                       massdensity_field ( k, joker, joker, joker ),
                       t_field, limits,
                       scat_data_meta_array, scat_data_start,
                       scat_data_nelem[k], part_species[k], delim,
                       verbosity);
    }

    //---- pnd_field calculations for MP48 -------------------------------
    else if ( psd_param == "MP48" )
    {
        psd = "MP48";

        //check for expected particle field name
        if ( partfield_name != "Rain" && partfield_name != "Snow")
        {
            out1 << "WARNING! The particle field name is unequal 'Rain' and 'Snow'.\n"
                 << psd << " should only be applied to liquid or frozen"
                 << " precipitation.\n";
        }

        //check for expected particle phase
        if ( part_type != "Ice" &&  part_type != "Water" )
        {
            out1 << "WARNING! The particle phase is unequal 'Ice' and 'Water'.\n"
                 << psd << " should only be applied to liquid water or water"
                 << " ice particles.\n";
        }

        pnd_fieldMP48( pnd_field,
                       massdensity_field ( k, joker, joker, joker ),
                       limits,
                       scat_data_meta_array, scat_data_start,
                       scat_data_nelem[k], part_species[k], delim,
                       verbosity);
    }

    // ---- pnd_field calculations for liquid water clouds ---------------------
    //      (parameters from Hess98, Stratus continental)
    else if ( psd_param == "H98_STCO" || psd_param == "STCO" || psd_param == "liquid" )
    {
        psd = "H98_STCO";

        //check for expected particle field name
        if ( partfield_name != "LWC")
        {
            out1 << "WARNING! The particle field name is unequal 'LWC'.\n"
                 << psd << " should should only be applied to liquid or frozen"
                 << " precipitation.\n";
        }

        //check for expected particle phase
        if ( part_type != "Water" )
        {
            out1 << "WARNING! The particle phase is unequal 'Water'.\n"
                 << psd << " should only be applied to liquid water"
                 << " particles.\n";
        }

        pnd_fieldH98 (pnd_field,
                       massdensity_field ( k, joker, joker, joker ),
                       limits,
                       scat_data_meta_array, scat_data_start,
                       scat_data_nelem[k], part_species[k], delim,
                       verbosity);
    }

    else
    {
        ostringstream os;
        os << "Size distribution function '" << psd_param << "' is unknown!\n";
        throw runtime_error( os.str() );
    }

    // alter starting index of current scattering data array to starting index
    // of next iteration step
    scat_data_start = scat_data_start + scat_data_nelem[k];

  }
}

/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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
  \file   m_atmosphere.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-16

  \brief  Workspace functions to set variables defining the atmosphere.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "check_input.h"
#include "matpackI.h"
#include "messages.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! AtmosphereSet1D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void AtmosphereSet1D(
        // WS Output:
              Index&    atmosphere_dim,
              Vector&   lat_grid,
              Vector&   lon_grid )
{
  out2 << "  Sets the atmospheric dimensionality to 1.\n";
  out3 << "    atmosphere_dim = 1\n";
  out3 << "    lat_grid is set to be an empty vector\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  atmosphere_dim = 1;
  lat_grid.resize(0);
  lon_grid.resize(0);
}



//! AtmosphereSet2D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void AtmosphereSet2D(
        // WS Output:
              Index&    atmosphere_dim,
              Vector&   lon_grid,
              Numeric&  lat_1d,
              Numeric&  meridian_angle_1d )
{
  out2 << "  Sets the atmospheric dimensionality to 2.\n";
  out3 << "    atmosphere_dim = 2\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  out3 << "    lat_1d = -999\n";
  out3 << "    meridian_angle_1d = -999\n";
  atmosphere_dim = 2;
  lon_grid.resize(0);
  lat_1d = -999;
  meridian_angle_1d = -999;
}



//! AtmosphereSet3D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void AtmosphereSet3D(
        // WS Output:
              Index&    atmosphere_dim,
              Numeric&  latitude_1d,
              Numeric&  meridian_angle_1d )
{
  out2 << "  Sets the atmospheric dimensionality to 2.\n";
  out3 << "    atmosphere_dim = 3\n";
  out3 << "    lat_1d = -999\n";
  out3 << "    meridian_angle_1d = -999\n";
  atmosphere_dim = 3;
  latitude_1d = -999;
  meridian_angle_1d = -999;
}



//! e_groundSet
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-29
*/
void e_groundSet(
        // WS Output:
              Tensor3&   e_ground,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Vector&    f_grid,
        // Control Parameters:
        const Numeric& value )
{
  // Check input (use dummy for *p_grid*).
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, Vector(2,2,-1), lat_grid, lon_grid );
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );

  Index nlat = lat_grid.nelem();
  Index nlon = lon_grid.nelem();
  if( atmosphere_dim < 2 )
    { nlat = 1; }
  if( atmosphere_dim < 3 )
    { nlon = 1; }

  out2 << "  e_ground = " << value << "\n";
  out3 << "            npages : " << f_grid.nelem() << "\n";
  out3 << "            nrows  : " << nlat << "\n";
  out3 << "            ncols  : " << nlon << "\n";

  e_ground.resize( f_grid.nelem(), nlat, nlon );
  e_ground = value;
}



//! GroundSetToBlackbody
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-27
*/
void GroundSetToBlackbody(
        // WS Output:
              Index&     blackbody_ground,
              Tensor3&   e_ground )
{
  out2 << "  blackbody_ground = 1\n";
  out2 << "  e_ground = 1\n";
  out3 << "            npages : 1\n";
  out3 << "            nrows  : 1\n";
  out3 << "            ncols  : 1\n";
  blackbody_ground = 1;
  e_ground.resize(1,1,1);
  e_ground(0,0,0) = 1;
}



//! r_geoidWGS84
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void r_geoidWGS84(
        // WS Output:
              Matrix&    r_geoid,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Numeric&   latitude_1d,
        const Numeric&   azimuth_angle_1d )
{
  // Check input (use dummy for *p_grid*).
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, Vector(2,2,-1), lat_grid, lon_grid );

  // Radii of the reference ellipsiod and the values squared, where double
  // is used to avoid numerical problems.
  const Numeric rq  = 6378.138e3, rp  = 6356.752e3;
  const double  rq2 = rq*rq,      rp2 = rp*rp;

  // For 1D the geoid is set to curvature radius for the point and direction
  // given by latitude_1d and azimuth_angle_1d.
  if( atmosphere_dim == 1 )
    {
      chk_if_in_range( "latitude_1d", latitude_1d, -90, 90 );
      chk_if_in_range( "azimuth_angle_1d", azimuth_angle_1d, -180, 180 );

      out2 << "  Sets r_geoid to the curvature radius of the WGS-84 "
           << "reference ellipsiod.\n";

      r_geoid.resize(1,1);

      // Cosine and sine of the latitude. The values are only used squared.
      double cv = cos( latitude_1d * DEG2RAD ); 
             cv = cv * cv; 
      double sv = sin( latitude_1d * DEG2RAD );
             sv = sv * sv; 

      // Calculate NS and EW radius
      Numeric rns = rq2*rp2/pow(rq2*cv+rp2*sv,1.5);
      Numeric rew = rq2/sqrt(rq2*cv+rp2*sv);

      // Calculate the curvature radius in the observation direction
      cv = cos( azimuth_angle_1d * DEG2RAD );
      sv = sin( azimuth_angle_1d * DEG2RAD );
      r_geoid(0,0) = 1/(cv*cv/rns+sv*sv/rew);
    }

  else
    {
      out2 << "  Sets r_geoid to the radius of the WGS-84 "
           << "reference ellipsiod.\n";

      // Number of latitudes
      const Index nlat = lat_grid.nelem();

      // Effective number of longitudes
      Index nlon;
      if( atmosphere_dim == 2 )
        nlon = 1;
      else
        nlon = lon_grid.nelem();

      r_geoid.resize( nlat, nlon );

      Numeric sv, cv, rv;
      for( Index lat=0; lat<nlat; lat++ )
        {
	  // Check that the latutide is inside [-90,90]
          if( ( lat_grid[lat] < -90 ) || ( lat_grid[lat] > 90 ) )
	    {
	      ostringstream os;
	      os << "The accepted range for latitudes in this function is "
                 << "[-90,90],\nbut a value of " << lat_grid[lat] 
		 << " was found.";
	      throw runtime_error( os.str() );
	    }	    

          // Cosine and sine of the latitude.
          cv = cos( lat_grid[lat] * DEG2RAD ); 
          sv = sin( lat_grid[lat] * DEG2RAD );
      
          // The radius of the ellipsiod
          rv = sqrt( rq2*rp2 / ( rq2*sv*sv + rp2*cv*cv ) );
      
          // Loop longitudes and fill *r_geoid*.
          for( Index lon=0; lon<nlon; lon++ )
	    r_geoid(lat,lon) = rv;
        }
    }
  out3 << "            nrows  : " << r_geoid.nrows() << "\n";
  out3 << "            ncols  : " << r_geoid.ncols() << "\n";
}




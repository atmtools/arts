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
#include "agenda_class.h"
#include "arts.h"
#include "check_input.h"
#include "physics_funcs.h"
#include "matpackI.h"
#include "messages.h"
#include "rte.h"
#include "special_interp.h"

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



//! GroundNoScatteringSingleEmissivity
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-09-22
*/
void GroundNoScatteringSingleEmissivity(
              Matrix&    ground_emission, 
              Matrix&    ground_los, 
              Tensor4&   ground_refl_coeffs,
        const Vector&    f_grid,
	const Index&     stokes_dim,
	const GridPos&   a_gp_p,
	const GridPos&   a_gp_lat,
	const GridPos&   a_gp_lon,
	const Vector&    a_los,
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
	const Matrix&    r_geoid,
        const Matrix&    z_ground,
        const Tensor3&   t_field,
        const Numeric&   e )
{
  //--- Check input -----------------------------------------------------------
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_ground", z_ground, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  chk_if_in_range( "the keyword *e*", e, 0, 1 );
  //---------------------------------------------------------------------------

  out2 << 
        "  Sets the ground to have no scattering and a constant emissivity.\n";
  out3 << 
        "     Ground temperature is obtained by interpolation of *t_field*.\n";
  out3 << "     Zenith angle for upwelling radiation is   : " 
       << a_los[0] << "\n";
  if( atmosphere_dim > 2 )
    {
      out3 << "     Azimuth angle for upwelling radiation is  : " << a_los[1] 
           << "\n";
    }

  // Some sizes
  const Index   nf = f_grid.nelem();

  // Resize output arguments
  ground_emission.resize( nf, stokes_dim );
  ground_los.resize( 1, a_los.nelem() );
  ground_refl_coeffs.resize( 1, nf, stokes_dim, stokes_dim );
  ground_refl_coeffs = 0;

  // Determine the temperature at the point of the ground reflection
  const Numeric t = interp_atmfield_by_gp( atmosphere_dim, p_grid, lat_grid, 
                    lon_grid, t_field, "t_field", a_gp_p, a_gp_lat, a_gp_lon );

  // Fill ground_emission and ground_refl_coeffs
  for( Index i=0; i<nf; i++ )
    { 
      ground_emission(i,0) = e * planck( f_grid[i], t ); 
      ground_refl_coeffs(0,i,0,0) = 1 - e;
      for( Index is=1; is<stokes_dim; is++ )
	{ 
	  ground_emission(i,is) = 0; 
	  ground_refl_coeffs(0,i,is,is) = 1 - e;
	}
      // Note that other elements of ground_feld_coeffs are set to 0 above.
    }

  // Calculate LOS for downwelling radiation
  ground_specular_los( ground_los(0,Range(joker)), atmosphere_dim, r_geoid,
                      z_ground, lat_grid, lon_grid, a_gp_lat,a_gp_lon, a_los );

  out3 << "     Zenith angle for downwelling radiation is : " 
       << ground_los(0,0) << "\n";
  if( atmosphere_dim > 2 )
    {
      out3 << "     Azimuth angle for downwelling radiation is: " 
           << ground_los(0,1) << "\n";
    }
}



//! GroundTreatAsBlackbody
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-09-17
*/
void GroundTreatAsBlackbody(
              Matrix&    ground_emission, 
              Matrix&    ground_los, 
              Tensor4&   ground_refl_coeffs,
        const Vector&    f_grid,
	const Index&     stokes_dim,
	const GridPos&   a_gp_p,
	const GridPos&   a_gp_lat,
	const GridPos&   a_gp_lon,
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Tensor3&   t_field )
{
  //--- Check input -----------------------------------------------------------
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  //---------------------------------------------------------------------------

  out2 << "  Sets the ground to be a blackbody.\n";
  out3 << 
        "     Ground temperature is obtained by interpolation of *t_field*.\n";

  // Set ground_los and ground_refl_coeffs to be empty as no downwelling
  // spectra shall be calculated
  ground_los.resize(0,0);
  ground_refl_coeffs.resize(0,0,0,0);
  
  // Some sizes
  const Index   nf = f_grid.nelem();

  // Resize ground_emission.
  ground_emission.resize(nf,stokes_dim);

  // Determine the temperature at the point of the ground reflection
  const Numeric t = interp_atmfield_by_gp( atmosphere_dim, p_grid, lat_grid, 
                    lon_grid, t_field, "t_field", a_gp_p, a_gp_lat, a_gp_lon );

  // Fill ground_emission with unpolarised blackbody radiation
  for( Index i=0; i<nf; i++ )
    { 
      ground_emission(i,0) = planck( f_grid[i], t ); 
      for( Index is=1; is<stokes_dim; is++ )
	{ ground_emission(i,is) = 0; }
    }
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




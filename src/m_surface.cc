/* Copyright (C) 2008
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
  \date   2008-09-17

  \brief  Workspace functions associated with the geoid and the surface.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/




/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
//#include <stdexcept>
#include "arts.h"
#include "check_input.h"
#include "complex.h"          
#include "physics_funcs.h"
#include "matpackIII.h"
#include "math_funcs.h"
#include "messages.h"
#include "special_interp.h"
#include "absorption.h"
#include "interpolation.h"
#include "fastem.h"
#include "gridded_fields.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric EARTH_RADIUS;

extern const Index GFIELD4_IA_GRID;
extern const Index GFIELD4_F_GRID;
extern const Index GFIELD4_LAT_GRID;
extern const Index GFIELD4_LON_GRID;




/*===========================================================================
  === Help functions (not WSMs)
  ===========================================================================*/

//! surface_specular_los
/*!
    Calculates the LOS for a specular surface reflection.

    \param   los                In/Out: LOS to be modified.
    \param   atmosphere_dim     Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-09-22
*/
void surface_specular_los(
              VectorView   los,
        const Index&       atmosphere_dim )
{
  assert( atmosphere_dim >= 1  &&  atmosphere_dim <= 3 );

  if( atmosphere_dim == 1 )
    {
      assert( los.nelem() == 1 );
      assert( los[0] > 90 );      // Otherwise surface refl. not possible
      assert( los[0] <= 180 ); 

      los[0] = 180 - los[0];
    }

  else if( atmosphere_dim == 2 )
    {
      assert( los.nelem() == 1 );
      assert( abs(los[0]) <= 180 ); 

      los[0] = sign( los[0] ) * 180 - los[0];
    }

  else if( atmosphere_dim == 3 )
    {
      assert( los.nelem() == 2 );
      assert( los[0] >= 0 ); 
      assert( los[0] <= 180 ); 
      assert( abs( los[1] ) <= 180 ); 

      // Calculate LOS neglecting any tilt of the surface
      los[0] = 180 - los[0];
      //los[1] = los[1];
    }
}



//! surface_specular_R_and_b
/*!
    Sets up the surface reflection matrix and emission vector for
    the case of specular reflection.

    The function handles only one frequency at the time.

    See further the surface chapter in the user guide.

    \param   surface_rmatrix   Out: As the WSV with the same name, but slice
                                    for one direction and one frequency.
    \param   surface_emission  Out: As the WSV with the same name, but slice
                                    for one direction and one frequency.
    \param   Rv                In: Complex amplitude relection coefficient
                                   for vertical polarisation.
    \param   Rh                In: Complex amplitude relection coefficient
                                   for horisontal polarisation.
    \param   f                 In: Frequency (a scalar).
    \param   stokes_dim        In: As the WSV with the same name.
    \param   surface_skin_t    In: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2004-09-24
*/
void surface_specular_R_and_b(
              MatrixView   surface_rmatrix,
              VectorView   surface_emission,
        const Complex&     Rv,
        const Complex&     Rh,
        const Numeric&     f,
        const Index&       stokes_dim,
        const Numeric&     surface_skin_t )
{
  assert( surface_rmatrix.nrows() == stokes_dim );
  assert( surface_rmatrix.ncols() == stokes_dim );
  assert( surface_emission.nelem() == stokes_dim );

  // Expressions are derived in the surface chapter in the user guide

  const Numeric   rv    = pow( abs(Rv), 2.0 );
  const Numeric   rh    = pow( abs(Rh), 2.0 );
  const Numeric   rmean = ( rv + rh ) / 2;
  const Numeric   B     = planck( f, surface_skin_t );

  surface_rmatrix   = 0.0;
  surface_emission  = 0.0;

  surface_rmatrix(0,0) = rmean;
  surface_emission[0]  = B * ( 1 - rmean );

  if( stokes_dim > 1 )
    {
      const Numeric   rdiff = ( rv - rh ) / 2;

      surface_rmatrix(1,0) = rdiff;
      surface_rmatrix(0,1) = rdiff;
      surface_rmatrix(1,1) = rmean;
      surface_emission[1]  = -B * rdiff;

        if( stokes_dim > 2 )
          {
            const Complex   a     = Rh * conj(Rv);
            const Complex   b     = Rv * conj(Rh);
            const Numeric   c     = real( a + b ) / 2.0;

            surface_rmatrix(2,2) = c;
      
            if( stokes_dim > 3 )
              {
                const Numeric   d     = imag( a - b ) / 2.0;

                surface_rmatrix(3,2) = d;
                surface_rmatrix(2,3) = d;
                surface_rmatrix(3,3) = c;
              }
          }
    }
}





/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void InterpSurfaceFieldToRteGps(
                 Numeric&   outvalue,
           const Index&     atmosphere_dim,
           const GridPos&   rte_gp_lat,
           const GridPos&   rte_gp_lon,
           const Matrix&    field )
{
  // Interpolate
  outvalue = interp_atmsurface_by_gp( atmosphere_dim, field, 
                                      rte_gp_lat, rte_gp_lon );

  out3 << "    Result = " << outvalue << "\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void InterpSurfaceEmissivityFieldIncLatLon(
                 Numeric&   outvalue,
           const Index&     atmosphere_dim,
           const Vector&    rte_pos,
           const Vector&    rte_los,
           const GriddedField3&   gfield )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  if( rte_pos.nelem() != atmosphere_dim )
    throw runtime_error( "Length of *rte_pos* must match *atmoshere_dim*." );
    
  // Interpolate
  //
  Numeric lat = 0, lon = 0;
  if( atmosphere_dim >= 2 )
    lat = rte_pos[1];
  if( atmosphere_dim == 3 )
    lon = rte_pos[2];
  //
  interp_gfield3( outvalue, gfield, atmosphere_dim, 180-abs(rte_los[0]), lat, 
                  lon, "Incidence angle", "Latitude", "Longitude" );

  out3 << "    Result = " << outvalue << "\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void r_geoidSpherical(
        // WS Output:
              Matrix&    r_geoid,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Numeric&   r )
{
  // Check input (use dummy for *p_grid*).
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, Vector(2,2,-1), lat_grid, lon_grid );

  // What radius to use?
  Numeric   r_local = r;
  if( r < 0 )
    { r_local = EARTH_RADIUS; }

  out2 << "  Sets r_geoid to a sphere with a constant radius of " 
       << r_local/1e3 << " km.\n";

  // Effective number of latitudes and longitudes
  Index nlat=1, nlon=1;
  if( atmosphere_dim >= 2 )
    { nlat = lat_grid.nelem(); }
  if( atmosphere_dim >= 3 )
    { nlon = lon_grid.nelem(); }
        
  r_geoid.resize( nlat, nlon );

  r_geoid = r_local;
  
  out3 << "            nrows  : " << r_geoid.nrows() << "\n";
  out3 << "            ncols  : " << r_geoid.ncols() << "\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void r_geoidWGS84(
        // WS Output:
              Matrix&    r_geoid,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid )
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
      const Numeric azimuth_angle_1d = 0.0;

      if( lat_grid.nelem() != 1 )
        throw runtime_error( 
               "This method requires, for 1D, that *lat_grid* has length 1." );

      chk_if_in_range( "element of lat_grid", lat_grid[0], -90., 90. );

      out2 << "  Sets r_geoid to the curvature radius of the WGS-84 "
           << "reference ellipsiod.\n";

      r_geoid.resize(1,1);

      // Cosine and sine of the latitude. The values are only used squared.
      double cv = cos( lat_grid[0] * DEG2RAD ); 
             cv = cv * cv; 
      double sv = sin( lat_grid[0] * DEG2RAD );
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



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceBlackbody(
              Matrix&    surface_los,
              Tensor4&   surface_rmatrix,
              Matrix&    surface_emission,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const Numeric&   surface_skin_t )
{
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );

  out2 << "  Sets variables to model a blackbody surface with a temperature "
       << " of " << surface_skin_t << " K.\n";
  surface_los.resize(0,0);
  surface_rmatrix.resize(0,0,0,0);

  const Index   nf = f_grid.nelem();
  surface_emission.resize( nf, stokes_dim );
  surface_emission = 0.0;
  for( Index iv=0; iv<nf; iv++ )
    { 
      surface_emission(iv,0) = planck( f_grid[iv], surface_skin_t );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatReflectivity(
              Matrix&    surface_los,
              Tensor4&   surface_rmatrix,
              Matrix&    surface_emission,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const Index&     atmosphere_dim,
        const Vector&    rte_los,
        const Numeric&   surface_skin_t,
        const Vector&    surface_scalar_reflectivity )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );

  const Index   nf = f_grid.nelem();

  if( surface_scalar_reflectivity.nelem() != nf  &&  
      surface_scalar_reflectivity.nelem() != 1 )
    {
      ostringstream os;
      os << "The number of elements in *surface_scalar_reflectivity* should\n"
         << "match length of *f_grid* or be 1."
         << "\n length of *f_grid* : " << nf 
         << "\n length of *surface_scalar_reflectivity* : " 
         << surface_scalar_reflectivity.nelem()
         << "\n";
      throw runtime_error( os.str() );
    }

  if( min(surface_scalar_reflectivity) < 0  ||  
      max(surface_scalar_reflectivity) > 1 )
    {
      throw runtime_error( 
         "All values in *surface_scalar_reflectivity* must be inside [0,1]." );
    }

  out2 << "  Sets variables to model a flat surface\n";
  out3 << "     surface temperature: " << surface_skin_t << " K.\n";

  surface_los.resize( 1, rte_los.nelem() );
  surface_los(0,joker) = rte_los;
  surface_specular_los( surface_los(0, joker) , atmosphere_dim );

  surface_emission.resize( nf, stokes_dim );
  surface_rmatrix.resize(1,nf,stokes_dim,stokes_dim);
  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  Numeric r = 0.0;

  for( Index iv=0; iv<nf; iv++ )
    { 
      if( iv == 0  || surface_scalar_reflectivity.nelem() > 1 )
        { r = surface_scalar_reflectivity[iv]; }

      surface_emission(iv,0) = (1.0-r) * planck( f_grid[iv], surface_skin_t );
      for( Index is=0; is<stokes_dim; is++ )
        { surface_rmatrix(0,iv,is,is) = r; }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatRefractiveIndex(
              Matrix&    surface_los,
              Tensor4&   surface_rmatrix,
              Matrix&    surface_emission,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const Index&     atmosphere_dim,
        const Vector&    rte_los,
        const Numeric&   surface_skin_t,
        const Matrix&    complex_n )
{
  const Index   nf = f_grid.nelem();

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );

  chk_matrix_ncols( "complex_n", complex_n, 2 ); 
  //
  if( complex_n.nrows() != nf  && complex_n.nrows() != 1 )
    {
      ostringstream os;
      os << "The number of rows in *complex_n* should be 1 or match the length "
         << "of *f_grid*,"
         << "\n length of *f_grid*  : " << nf 
         << "\n rows in *complex_n* : " << complex_n.nrows() << "\n";
      throw runtime_error( os.str() );
    }

  out2 << "  Sets variables to model a flat surface\n";
  out3 << "     surface temperature: " << surface_skin_t << " K.\n";

  surface_los.resize( 1, rte_los.nelem() );
  surface_los(0,joker) = rte_los;
  surface_specular_los( surface_los(0,joker), atmosphere_dim );

  surface_emission.resize( nf, stokes_dim );
  surface_rmatrix.resize( 1, nf, stokes_dim, stokes_dim );

  // Complex (amplitude) reflection coefficients
  Complex  Rv, Rh;

  for( Index iv=0; iv<nf; iv++ )
    { 
      if( iv == 0  || complex_n.nrows() > 1 )
        {
          // Set n2 (refractive index of surface medium)
          Complex n2( complex_n(iv,0), complex_n(iv,1) );

          // Amplitude reflection coefficients
          fresnel( Rv, Rh, Numeric(1.0), n2, 180.0-abs(rte_los[0]) );
        }

      // Fill reflection matrix and emission vector
      surface_specular_R_and_b( surface_rmatrix(0,iv,joker,joker), 
                                surface_emission(iv,joker), Rv, Rh, 
                                f_grid[iv], stokes_dim, surface_skin_t );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatVaryingEmissivity(
              Matrix&    surface_los,
              Tensor4&   surface_rmatrix,
              Matrix&    surface_emission,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const Index&     atmosphere_dim,
        const Vector&    rte_los,
        const Numeric&   surface_skin_t,
        const Vector&    surface_emissivity )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );

  const Index   nf = f_grid.nelem();

  if( surface_emissivity.nelem() != nf )
    {
      ostringstream os;
      os << "The number of elements in *surface_emissivity* should match\n"
         << "length of *f_grid*."
         << "\n length of *f_grid* : " << nf 
         << "\n length of *surface_emissivity* : " << surface_emissivity.nelem()
         << "\n";
      throw runtime_error( os.str() );
    }

  if( min(surface_emissivity) < 0  ||  max(surface_emissivity) > 1 )
    {
      throw runtime_error( 
                  "All values in *surface_emissivity* must be inside [0,1]." );
    }

  out2 << "  Sets variables to model a flat surface\n";
  out3 << "     surface temperature: " << surface_skin_t << " K.\n";

  surface_los.resize( 1, rte_los.nelem() );
  surface_los(0,joker) = rte_los;
  surface_specular_los( surface_los(0, joker) , atmosphere_dim );

  surface_emission.resize( nf, stokes_dim );
  surface_rmatrix.resize(1,nf,stokes_dim,stokes_dim);
  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  for( Index iv=0; iv<nf; iv++ )
    { 
      surface_emission(iv,0) = surface_emissivity[iv] * 
                                          planck( f_grid[iv], surface_skin_t );
      for( Index is=0; is<stokes_dim; is++ )
        { surface_rmatrix(0,iv,is,is) = 1 - surface_emissivity[iv]; }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatSingleEmissivity(
              Matrix&    surface_los,
              Tensor4&   surface_rmatrix,
              Matrix&    surface_emission,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const Index&     atmosphere_dim,
        const Vector&    rte_los,
        const Numeric&   surface_skin_t,
        const Numeric&   surface_emissivity )
{
  chk_if_in_range( "surface_emissivity", surface_emissivity, 0., 1. );

  Vector a_vector( f_grid.nelem(), surface_emissivity );

  surfaceFlatVaryingEmissivity( surface_los, surface_rmatrix, surface_emission, 
                                f_grid, stokes_dim, atmosphere_dim, rte_los, 
                                surface_skin_t, a_vector );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceLambertianSimple(
              Matrix&    surface_los,
              Tensor4&   surface_rmatrix,
              Matrix&    surface_emission,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const Index&     atmosphere_dim,
        const Vector&    rte_los,
        const Numeric&   surface_skin_t,
        const Vector&    surface_scalar_reflectivity,
        const Index&     np,
        const Numeric&   za_pos )
{
  const Index   nf = f_grid.nelem();

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 1 );
  chk_not_negative( "surface_skin_t", surface_skin_t );
  chk_if_in_range( "za_pos", za_pos, 0, 1 );

  if( surface_scalar_reflectivity.nelem() != 1  && 
      surface_scalar_reflectivity.nelem() != nf )
    {
      throw runtime_error( "Length of *surface_scalar_reflectivity* must "
                           "either be 1 or be equal to length of *f_grid*." );
    } 

  // Allocate and init everything to zero
  //
  surface_los.resize( np, rte_los.nelem() );
  surface_rmatrix.resize(np,nf,stokes_dim,stokes_dim);
  surface_emission.resize( nf, stokes_dim );
  //
  surface_los      = 0.0;
  surface_rmatrix  = 0.0;
  surface_emission = 0.0;

  // Help variables
  //
  const Numeric dza = 90.0 / np;
  const Vector za_lims( 0.0, np+1, dza );

  // surface_los
  for( Index ip=0; ip<np; ip++ )
    { surface_los(ip,0) = za_lims[ip] + za_pos * dza; }

  // Loop frequencies and set remaining values
  //
  Numeric rd = surface_scalar_reflectivity[0];
  chk_if_in_range( "surface_scalar_reflectivity[0]", rd, 0, 1 );
  //
  for( Index iv=0; iv<nf; iv++ )
    {
      // Update reflectivity?
      if( iv  &&  surface_scalar_reflectivity.nelem() > 1 )
        { 
          rd = surface_scalar_reflectivity[iv];
          chk_if_in_range( "surface_scalar_reflectivity[i]", rd, 0, 1 );
        }

      // surface_rmatrix
      for( Index ip=0; ip<np; ip++ )
        {
          const Numeric w = rd * 0.5 * ( cos(2*DEG2RAD*za_lims[ip]) - 
                                         cos(2*DEG2RAD*za_lims[ip+1]) );
          for( Index is=0; is<stokes_dim; is++ )
            { surface_rmatrix(ip,iv,is,is) = w; }
        }

      // surface_emission
      surface_emission(iv,0) = (1-rd) * planck( f_grid[iv], surface_skin_t );
    }
}

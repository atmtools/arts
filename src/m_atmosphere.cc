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
#include "matpackIII.h"
#include "messages.h"
#include "rte.h"
#include "special_interp.h"
#include "absorption.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "logic.h"
#include "xml_io.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric EARTH_RADIUS;



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



//! Calculate atmospheric fields.
/*!
  This method interpolates the data for atmospheric fields on the atmospheric
  grids used for the calculation.

   See also the the online help (arts -d FUNCTION_NAME)

   \author Claudia Emde
   \date   2002-11-29
   \date   2004-02-23 Modified. Used GriddedField3 instead of ArrayOfTensor3.
*/
void AtmFieldsCalc(//WS Output:
                   Tensor3& t_field,
                   Tensor3& z_field,
                   Tensor4& vmr_field,
                   //WS Input
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const GriddedField3& t_field_raw,
                   const GriddedField3& z_field_raw,
                   const ArrayOfGriddedField3& vmr_field_raw,
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
      t_field.resize(p_grid.nelem(), 1, 1);
      z_field.resize(p_grid.nelem(), 1, 1);
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), 1, 1);

 
      // Gridpositions:
      ArrayOfGridPos gp_p(p_grid.nelem());
  
      // Interpolate t_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, t_field_raw.p_grid, p_grid );

      // Interpolation weights:
      Matrix itw(p_grid.nelem(), 2);
      // (2 interpolation weights are required for 1D interpolation)
      interpweights( itw, gp_p);
  
      // Interpolate:
      interp( t_field(joker, 0, 0), itw, 
              t_field_raw.data(joker, 0, 0),  gp_p);

  
      // Interpolate z_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, z_field_raw.p_grid, p_grid );
     
      // Interpolation weights:
      interpweights( itw, gp_p );
      
      // Interpolate:
      interp( z_field(joker, 0, 0), itw,
              z_field_raw.data(joker, 0, 0), gp_p);
      
  
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, vmr_field_raw[gas_i].p_grid, p_grid);
      
          // Interpolation weights:
          interpweights( itw, gp_p);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, 0, 0),
                  itw, vmr_field_raw[gas_i].data(joker, 0, 0), gp_p);
        }
      
    }

  //=========================================================================
  else if(atmosphere_dim == 2)
    {
      if( t_field_raw.lat_grid.nelem() == 0 &&
          t_field_raw.lon_grid.nelem() == 0 )
        throw runtime_error(
                            "Raw data has wrong dimension (1D). "
                            "You have to use \n"
                            "AtmFieldsCalcExpand1D instead of AtmFieldsCalc."
                            );
      
      //Resize variables
      t_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
      z_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       1);
      
      
      // Gridpositions:
      ArrayOfGridPos gp_p(p_grid.nelem());
      ArrayOfGridPos gp_lat(lat_grid.nelem());
      
      // Interpolate t_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, t_field_raw.p_grid, p_grid );
      gridpos( gp_lat, t_field_raw.lat_grid, lat_grid );
            
      // Interpolation weights:
      Tensor3 itw(p_grid.nelem(), lat_grid.nelem(), 4);
      // (8 interpolation weights are required for 3D interpolation)
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( t_field(joker, joker, 0 ), itw,
              t_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
      
      
      // Interpolate z_field:
      // Calculate grid positions:
      p2gridpos( gp_p, z_field_raw.p_grid, p_grid );
      gridpos( gp_lat, z_field_raw.lat_grid, lat_grid );
            
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( z_field(joker, joker, 0), itw, 
              z_field_raw.data(joker, joker, 0), gp_p, gp_lat);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, vmr_field_raw[gas_i].p_grid, p_grid);
          gridpos(gp_lat, vmr_field_raw[gas_i].lat_grid, lat_grid);
                  
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, 0),
                  itw, vmr_field_raw[gas_i].data(joker, joker, 0),
                  gp_p, gp_lat);
        }
    }

  //================================================================
  // atmosphere_dim = 3    
  else
    {
      if( t_field_raw.lat_grid.nelem() == 0 &&
          t_field_raw.lon_grid.nelem() == 0 )
        throw runtime_error(
                            "Raw data has wrong dimension. You have to use \n"
                            "AtmFieldsCalcExpand1D instead of AtmFieldsCalc."
                            );

      //Resize variables
      t_field.resize(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
      z_field.resize(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       lon_grid.nelem());
      
      
      // Gridpositions:
      ArrayOfGridPos gp_p(p_grid.nelem());
      ArrayOfGridPos gp_lat(lat_grid.nelem());
      ArrayOfGridPos gp_lon(lon_grid.nelem());
      
      
      // Interpolate t_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, t_field_raw.p_grid, p_grid );
      gridpos( gp_lat, t_field_raw.lat_grid, lat_grid );
      gridpos( gp_lon, t_field_raw.lon_grid, lon_grid );
      
      // Interpolation weights:
      Tensor4 itw(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem(), 8);
      // (8 interpolation weights are required for 3D interpolation)
      interpweights( itw, gp_p, gp_lat, gp_lon );
      
      // Interpolate:
      interp( t_field, itw, t_field_raw.data,  gp_p, gp_lat, gp_lon);
      
      
      // Interpolate z_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, z_field_raw.p_grid, p_grid );
      gridpos( gp_lat, z_field_raw.lat_grid, lat_grid );
      gridpos( gp_lon, z_field_raw.lon_grid, lon_grid );
      
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat, gp_lon );
      
      // Interpolate:
      interp( z_field, itw, z_field_raw.data, gp_p, gp_lat, gp_lon);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, vmr_field_raw[gas_i].p_grid, p_grid);
          gridpos(gp_lat, vmr_field_raw[gas_i].lat_grid, lat_grid);
          gridpos(gp_lon, vmr_field_raw[gas_i].lon_grid, lon_grid);
          
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, joker),
                  itw, vmr_field_raw[gas_i].data, gp_p, gp_lat, gp_lon);
        }
    }
}




//! AtmFieldsCalcExpand1D
/*!

   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2003-01-09

   \date   2004-02-23 Modified by Claudia Emde:
           Used GriddedField3 instead of ArrayOfTensor3.
*/
void AtmFieldsCalcExpand1D(
            Tensor3&                 t_field,
            Tensor3&                 z_field,
            Tensor4&                 vmr_field,
      const Vector&                  p_grid,
      const Vector&                  lat_grid,
      const Vector&                  lon_grid,
      const GriddedField3&          t_field_raw,
      const GriddedField3&          z_field_raw,
      const ArrayOfGriddedField3&   vmr_field_raw,
      const Index&                   atmosphere_dim )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  if( atmosphere_dim == 1 )
    { throw runtime_error( 
     "This function is intended for 2D and 3D. For 1D, use *AtmFieldsCalc*.");}

  // Make 1D interpolation using some dummy variables
  Vector    vempty(0);
  Tensor3   t_temp, z_temp;
  Tensor4   vmr_temp;
  AtmFieldsCalc( t_temp, z_temp, vmr_temp, p_grid, vempty, vempty, 
                                  t_field_raw, z_field_raw, vmr_field_raw, 1 );

  // Move values from the temporary tensors to the return arguments
  const Index   np = p_grid.nelem();
  const Index   nlat = lat_grid.nelem();
        Index   nlon = lon_grid.nelem();
  if( atmosphere_dim == 2 )
    { nlon = 1; }
  const Index   nspecies = vmr_temp.nbooks();
  //
  assert( t_temp.npages() == np );
  //
  t_field.resize( np, nlat, nlon );
  z_field.resize( np, nlat, nlon );
  vmr_field.resize( nspecies, np, nlat, nlon );
  //
  for( Index ilon=0; ilon<nlon; ilon++ )
    {
      for( Index ilat=0; ilat<nlat; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              t_field(ip,ilat,ilon) = t_temp(ip,0,0);
              z_field(ip,ilat,ilon) = z_temp(ip,0,0);
              for( Index is=0; is<nspecies; is++ )
                { vmr_field(is,ip,ilat,ilon) = vmr_temp(is,ip,0,0); }
            }
        }
    }
}


  
//! Read atmospheric scenario.
/* 

   See the the online help (arts -d FUNCTION_NAME)

 \param t_field_raw temperature field data
 \param z_field_raw altitude field data
 \param vmr_field_raw vmr field data
 \param gas_species gas species for calculation
 \param basename name of scenario

 \author Claudia Emde
 \date   2002-11-29
*/
void AtmRawRead(//WS Output:
                GriddedField3& t_field_raw,
                GriddedField3& z_field_raw,
                ArrayOfGriddedField3& vmr_field_raw,
                //WS Input:
                const ArrayOfArrayOfSpeciesTag& gas_species,
                //Keyword:
                const String& basename)
{
  // Read the temperature field:
  String file_name = basename + ".t.xml";
  xml_read_from_file( file_name, t_field_raw);
  
  out3 << "Temperature field read from file: " << file_name << "\n";  

  // Read geometrical altitude field:
  file_name = basename + ".z.xml";
  xml_read_from_file( file_name, z_field_raw);

  out3 << "Altitude field read from file: " << file_name << "\n";  


  // The species lookup data:

  extern const Array<SpeciesRecord> species_data;
  
  // We need to read one profile for each tag group.
  for ( Index i=0; i<gas_species.nelem(); i ++)
    {
      // Determine the name.
      file_name =
        basename + "." +
        species_data[gas_species[i][0].Species()].Name() + ".xml";
      
      // Add an element for this tag group to the vmr profiles:
      GriddedField3 vmr_field_data;
      vmr_field_raw.push_back(vmr_field_data);
      
      // Read the VMR:
      xml_read_from_file( file_name, vmr_field_raw[vmr_field_raw.nelem()-1]);
      
      // state the source of profile.
      out3 << "  " << species_data[gas_species[i][0].Species()].Name()
           << " profile read from file: " << file_name << "\n";
    }
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
  out2 << "  Sets the atmospheric dimensionality to 3.\n";
  out3 << "    atmosphere_dim = 3\n";
  out3 << "    lat_1d = -999\n";
  out3 << "    meridian_angle_1d = -999\n";
  atmosphere_dim = 3;
  latitude_1d = -999;
  meridian_angle_1d = -999;
}



//! surfaceNoScatteringSingleEmissivity
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-09-22
*/
void surfaceNoScatteringSingleEmissivity(
              Matrix&    surface_emission, 
              Matrix&    surface_los, 
              Tensor4&   surface_refl_coeffs,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const GridPos&   rte_gp_p,
        const GridPos&   rte_gp_lat,
        const GridPos&   rte_gp_lon,
        const Vector&    rte_los,
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Matrix&    r_geoid,
        const Matrix&    z_surface,
        const Tensor3&   t_field,
        const Numeric&   e )
{
  //--- Check input -----------------------------------------------------------
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_surface", z_surface, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  chk_if_in_range( "the keyword *e*", e, 0, 1 );
  //---------------------------------------------------------------------------

  out2 << 
        "  Sets the surface to have no scattering and a constant emissivity.\n";
  out3 << 
        "     Surface temperature is obtained by interpolation of *t_field*.\n";

  // Some sizes
  const Index   nf = f_grid.nelem();

  // Resize output arguments
  surface_emission.resize( nf, stokes_dim );
  surface_los.resize( 1, rte_los.nelem() );
  surface_refl_coeffs.resize( 1, nf, stokes_dim, stokes_dim );
  surface_refl_coeffs = 0;

  // Determine the temperature at the point of the surface reflection
  const Numeric t = interp_atmfield_by_gp( atmosphere_dim, p_grid, lat_grid, 
              lon_grid, t_field, "t_field", rte_gp_p, rte_gp_lat, rte_gp_lon );

  out3 << "     Surface temperature is                     : " << t << "\n";

  // Fill surface_emission and surface_refl_coeffs
  for( Index i=0; i<nf; i++ )
    { 
      surface_emission(i,0) = e * planck( f_grid[i], t ); 
      surface_refl_coeffs(0,i,0,0) = 1 - e;
      for( Index is=1; is<stokes_dim; is++ )
        { 
          surface_emission(i,is) = 0; 
          surface_refl_coeffs(0,i,is,is) = 1 - e;
        }
      // Note that other elements of surface_refl_coeffs are set to 0 above.
    }

  // Calculate LOS for downwelling radiation
  surface_specular_los( surface_los(0,joker), atmosphere_dim, r_geoid,
                z_surface, lat_grid, lon_grid, rte_gp_lat,rte_gp_lon, rte_los );

  out3 << "     Zenith angle for upwelling radiation is   : " 
       << rte_los[0] << "\n";
  if( atmosphere_dim > 2 )
    {
      out3 << "     Azimuth angle for upwelling radiation is  : " 
           << rte_los[1] << "\n";
    }
  out3 << "     Zenith angle for downwelling radiation is : " 
       << surface_los(0,0) << "\n";
  if( atmosphere_dim > 2 )
    {
      out3 << "     Azimuth angle for downwelling radiation is: " 
           << surface_los(0,1) << "\n";
    }
}



//! surfaceTreatAsBlackbody
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-09-17
*/
void surfaceTreatAsBlackbody(
              Matrix&    surface_emission, 
              Matrix&    surface_los, 
              Tensor4&   surface_refl_coeffs,
        const Vector&    f_grid,
        const Index&     stokes_dim,
        const GridPos&   rte_gp_p,
        const GridPos&   rte_gp_lat,
        const GridPos&   rte_gp_lon,
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

  out2 << "  Sets the surface to be a blackbody.\n";
  out3 << 
        "     Surface temperature is obtained by interpolation of *t_field*.\n";

  // Set surface_los and surface_refl_coeffs to be empty as no downwelling
  // spectra shall be calculated
  surface_los.resize(0,0);
  surface_refl_coeffs.resize(0,0,0,0);
  
  // Some sizes
  const Index   nf = f_grid.nelem();

  // Resize surface_emission.
  surface_emission.resize(nf,stokes_dim);

  // Determine the temperature at the point of the surface reflection
  const Numeric t = interp_atmfield_by_gp( atmosphere_dim, p_grid, lat_grid, 
              lon_grid, t_field, "t_field", rte_gp_p, rte_gp_lat, rte_gp_lon );

  // Fill surface_emission with unpolarised blackbody radiation
  for( Index i=0; i<nf; i++ )
    { 
      surface_emission(i,0) = planck( f_grid[i], t ); 
      for( Index is=1; is<stokes_dim; is++ )
        { surface_emission(i,is) = 0; }
    }
}



//! r_geoidSpherical
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-12_02
*/
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




/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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
  \file   m_surface.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2008-09-17

  \brief  Workspace functions associated wih the surface and its properties.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/




/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "complex.h"          
#include "geodetic.h"          
#include "math_funcs.h"
#include "messages.h"
#include "special_interp.h"
#include "interpolation.h"
#include "physics_funcs.h"
#include "rte.h"
#include "surface.h"

extern const Numeric DEG2RAD;




/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void complex_nFromGriddedField4(
          Matrix&          complex_n,
    const Index&           stokes_dim,
    const Vector&          f_grid,
    const Index&           atmosphere_dim,
    const Vector&          lat_grid,
    const Vector&          lat_true,
    const Vector&          lon_true,
    const Vector&          rte_pos,
    const Vector&          rte_los,
    const GriddedField4&   n_field,
    const Verbosity&)
{
  // Set expected order of grids
  Index gfield_fID = 0;
  Index gfield_compID = 1;
  Index gfield_latID = 2;
  Index gfield_lonID = 3;

  // Basic checks and sizes
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_latlon_true( atmosphere_dim, lat_grid, lat_true, lon_true );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_los( atmosphere_dim, rte_los );
  n_field.checksize_strict();
  //
  chk_griddedfield_gridname( n_field, gfield_fID, "Frequency" );
  chk_griddedfield_gridname( n_field, gfield_compID, "Complex" );
  chk_griddedfield_gridname( n_field, gfield_latID, "Latitude" );
  chk_griddedfield_gridname( n_field, gfield_lonID, "Longitude" );
  //
  const Index nf_in = n_field.data.nbooks();
  const Index nn    = n_field.data.npages();
  const Index nlat  = n_field.data.nrows();
  const Index nlon  = n_field.data.ncols();
  //
  if( nlat < 2  ||  nlon < 2 )
    {
      ostringstream os;
      os << "The data in *n_field* must span a geographical region. That is,\n"
         << "the latitude and longitude grids must have a length >= 2.";
    } 
  //
  if( nn != 2 )
    {
      ostringstream os;
      os << "The data in *n_field* must have exactly two pages. One page each\n"
         << "for the real and the imaginary part of the complex refractive index.";
    } 

  const Vector& GFlat=n_field.get_numeric_grid(gfield_latID);
  const Vector& GFlon=n_field.get_numeric_grid(gfield_lonID);

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon( lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, 
                                                           lon_true, rte_pos );

  // Ensure correct coverage of lon grid
  Vector lon_shifted;
  lon_shiftgrid( lon_shifted, GFlon, lon[0] );

  // Check if lat/lon we need are actually covered
  chk_if_in_range( "rte_pos.lat", lat[0], GFlat[0],
                    GFlat[nlat-1] );
  chk_if_in_range( "rte_pos.lon", lon[0], lon_shifted[0],
                    lon_shifted[nlon-1] );

  // Interpolate in lat and lon
  Matrix n_f( nf_in, nn );
  {
    GridPos gp_lat, gp_lon;
    gridpos( gp_lat, GFlat, lat[0] );
    gridpos( gp_lon, lon_shifted, lon[0] );
    Vector itw(4);
    interpweights( itw, gp_lat, gp_lon );
    for( Index iv=0; iv<nf_in; iv++ )
      {
        for( Index in=0; in<nn; in++ )
          { 
            n_f(iv,in) = interp( itw, n_field.data(iv,in,joker,joker), 
                                                              gp_lat, gp_lon );
          }
      } 

  }    
  
  // Extract or interpolate in frequency
  //
  if( nf_in == 1 )
    { complex_n = n_f; }
  else
    {
      const Index nf_out = f_grid.nelem();
      complex_n.resize( nf_out, 2 );
      //
      ArrayOfGridPos gp( nf_out );
      Matrix         itw( nf_out, 2 );
      gridpos( gp, n_field.get_numeric_grid(gfield_fID), f_grid );
      interpweights( itw, gp );
      interp( complex_n(joker,0), itw, n_f(joker,0), gp );
      interp( complex_n(joker,1), itw, n_f(joker,1), gp );
    }     
}



/* Workspace method: Doxygen documentation will be auto-generated */
void InterpSurfaceFieldToRtePos(
          Numeric&   outvalue,
    const Index&     atmosphere_dim,
    const Vector&    lat_grid,
    const Vector&    lon_grid,
    const Vector&    pos,
    const Matrix&    field,
    const Verbosity& verbosity)
{
  // Input checks (dummy p_grid)
  chk_atm_grids( atmosphere_dim, Vector(2,2,-1), lat_grid, lon_grid );
  chk_atm_surface( "input argument *field*", field, atmosphere_dim, lat_grid, 
                                                                    lon_grid );
  chk_rte_pos( atmosphere_dim, pos );

  if( atmosphere_dim == 1 )
    { outvalue = field(0,0); }
  else
    {      
      chk_interpolation_grids( "Latitude interpolation", lat_grid, pos[1] );
      GridPos gp_lat, gp_lon;
      gridpos( gp_lat, lat_grid, pos[1] );
      if( atmosphere_dim == 3 )
        { 
          chk_interpolation_grids( "Longitude interpolation", lon_grid, pos[2]);
          gridpos( gp_lon, lon_grid, pos[2] );
        }
      //
      outvalue = interp_atmsurface_by_gp( atmosphere_dim, field, gp_lat, 
                                                                 gp_lon );
    }
  
  // Interpolate
  CREATE_OUT3;
  out3 << "    Result = " << outvalue << "\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceRtpropAgenda(
          Workspace&        ws,
          Matrix&           iy,
          Matrix&           iy_error,
          Index&            iy_error_type,
          Matrix&           iy_aux,
          ArrayOfTensor3&   diy_dx,  
    const Tensor3&          iy_transmission,
    const Index&            jacobian_do,
    const Index&            atmosphere_dim,
    const Vector&           lat_grid,
    const Vector&           lon_grid,
    const Vector&           refellipsoid,
    const Matrix&           z_surface,          
    const Tensor3&          t_field,
    const Tensor3&          z_field,
    const Tensor4&          vmr_field,
    const Index&            cloudbox_on,
    const Index&            stokes_dim,
    const Vector&           f_grid,
    const Vector&           rte_pos,
    const Vector&           rte_los,
    const Agenda&           iy_clearsky_agenda,
    const Agenda&           surface_rtprop_agenda,
    const Verbosity& )
{
  // Input checks
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_los( atmosphere_dim, rte_los );

  // Calculate specular direction
  Vector specular_los( max(Index(1),atmosphere_dim-1) );
  if( atmosphere_dim == 1 )
    { 
      chk_if_in_range( "1D rte_los zenith angle", rte_los[0], 90, 180 );
      specular_los[0] = 180 - rte_los[0];
    }
  else if( atmosphere_dim == 2 )
    { 
      chk_interpolation_grids( "Latitude interpolation", lat_grid, rte_pos[1] );
      GridPos gp_lat;
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      Numeric c1;         // Radial slope of the surface
      plevel_slope_2d( c1, lat_grid, refellipsoid, z_surface(joker,0), 
                                                          gp_lat, rte_los[0] );
      Vector itw(2); interpweights( itw, gp_lat );
      const Numeric rv_surface = refell2d( refellipsoid, lat_grid, gp_lat )
                                + interp( itw, z_surface(joker,0), gp_lat );
      specular_los[0] = 180-rte_los[0]-2*plevel_angletilt( rv_surface, c1 );
    }
  else if( atmosphere_dim == 3 )
    { 
      // Calculate surface normal in South-North direction
      chk_interpolation_grids( "Latitude interpolation", lat_grid, rte_pos[1] );
      chk_interpolation_grids( "Longitude interpolation", lon_grid, rte_pos[2]);
      GridPos gp_lat, gp_lon;
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      gridpos( gp_lon, lon_grid, rte_pos[2] );
      Numeric c1, c2;
      plevel_slope_3d( c1, c2, lat_grid, lon_grid, refellipsoid, z_surface, 
                       gp_lat, gp_lon, 0 );
      Vector itw(4); interpweights( itw, gp_lat, gp_lon );
      const Numeric rv_surface = refell2d( refellipsoid, lat_grid, gp_lat )
                                 + interp( itw, z_surface, gp_lat, gp_lon );
      const Numeric zaSN = 90 - plevel_angletilt( rv_surface, c1 );
      // The same for East-West
      plevel_slope_3d( c1, c2, lat_grid, lon_grid, refellipsoid, z_surface, 
                       gp_lat, gp_lon, 90 );
      const Numeric zaEW = 90 - plevel_angletilt( rv_surface, c1 );
      // Convert to Cartesian, and determine normal by cross-product
      Vector tangentSN(3), tangentEW(3), normal(3);
      zaaa2cart( tangentSN[0], tangentSN[1], tangentSN[2], zaSN, 0 );
      zaaa2cart( tangentEW[0], tangentEW[1], tangentEW[2], zaEW, 90 );
      cross3( normal, tangentSN, tangentEW );
   }


  // Call *surface_rtprop_agenda*
  Matrix    surface_los;
  Tensor4   surface_rmatrix;
  Matrix    surface_emission;
  //
  surface_rtprop_agendaExecute( ws, surface_emission, surface_los, 
                                surface_rmatrix, rte_pos, rte_los, 
                                surface_rtprop_agenda );

  // Check output of *surface_rtprop_agenda*
  const Index   nlos = surface_los.nrows();
  const Index   nf   = f_grid.nelem();
  //
  if( nlos )   // if 0, blackbody ground and not all checks are needed
    {
      if( surface_los.ncols() != rte_los.nelem() )
        throw runtime_error( 
                        "Number of columns in *surface_los* is not correct." );
      if( nlos != surface_rmatrix.nbooks() )
        throw runtime_error( 
                  "Mismatch in size of *surface_los* and *surface_rmatrix*." );
      if( surface_rmatrix.npages() != nf )
        throw runtime_error( 
                       "Mismatch in size of *surface_rmatrix* and *f_grid*." );
      if( surface_rmatrix.nrows() != stokes_dim  ||  
          surface_rmatrix.ncols() != stokes_dim ) 
        throw runtime_error( 
              "Mismatch between size of *surface_rmatrix* and *stokes_dim*." );
    }
  if( surface_emission.ncols() != stokes_dim )
    throw runtime_error( 
             "Mismatch between size of *surface_emission* and *stokes_dim*." );
  if( surface_emission.nrows() != nf )
    throw runtime_error( 
                       "Mismatch in size of *surface_emission* and f_grid*." );

  // Variable to hold down-welling radiation
  Tensor3   I( nlos, nf, stokes_dim );
 
  // Loop *surface_los*-es. If no such LOS, we are ready.
  if( nlos > 0 )
    {
      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          Vector los = surface_los(ilos,joker);

          // Include surface reflection matrix in *iy_transmission*
          // If iy_transmission is empty, this is interpreted as the
          // variable is not needed.
          //
          Tensor3 iy_trans_new;
          //
          if( iy_transmission.npages() )
            {
              iy_transmission_mult(  iy_trans_new, iy_transmission, 
                                     surface_rmatrix(ilos,joker,joker,joker) );
            }

          // Calculate downwelling radiation for LOS ilos 
          //
          // The variable iy_clearsky_agenda can in fact be 
          // iy_clearsky_BASIC_agenda
          //
          if( iy_clearsky_agenda.name() == "iy_clearsky_basic_agenda" )
            {
              iy_clearsky_basic_agendaExecute( ws, iy, rte_pos, los,
                                              cloudbox_on, iy_clearsky_agenda);
            }
          else
            {
              iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                                  iy_aux, diy_dx, 0, iy_trans_new, rte_pos, 
                                  los, cloudbox_on, jacobian_do, t_field,
                                  z_field, vmr_field, -1, iy_clearsky_agenda );
            }

          if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
            {
              ostringstream os;
              os << "The size of *iy* returned from *" 
                 << iy_clearsky_agenda.name() << "* is\n"
                 << "not correct:\n"
                 << "  expected size = [" << nf << "," << stokes_dim << "]\n"
                 << "  size of iy    = [" << iy.nrows() << "," << iy.ncols()
                 << "]\n";
              throw runtime_error( os.str() );      
            }

          I(ilos,joker,joker) = iy;
        }
    }
  
  // Add up
  surface_calc( iy, I, surface_los, surface_rmatrix, surface_emission );
}




/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceBlackbody(
          Workspace& ws,
          Matrix&    surface_los,
          Tensor4&   surface_rmatrix,
          Matrix&    surface_emission,
    const Vector&    f_grid,
    const Index&     stokes_dim,
    const Numeric&   surface_skin_t,
    const Agenda&    blackbody_radiation_agenda,
    const Verbosity& verbosity)
{  
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );

  CREATE_OUT2;
  out2 << "  Sets variables to model a blackbody surface with a temperature "
       << " of " << surface_skin_t << " K.\n";

  surface_los.resize(0,0);
  surface_rmatrix.resize(0,0,0,0);

  Vector b;
  blackbody_radiation_agendaExecute( ws, b, surface_skin_t, f_grid, 
                                     blackbody_radiation_agenda ); 

  const Index   nf = f_grid.nelem();
  surface_emission.resize( nf, stokes_dim );
  surface_emission = 0.0;

  for( Index iv=0; iv<nf; iv++ )
    { 
      surface_emission(iv,0) = b[iv]; 
      for( Index is=1; is<stokes_dim; is++ )
        { surface_emission(iv,is) = 0; } 
   }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatRefractiveIndex(
          Workspace& ws,
          Matrix&    surface_los,
          Tensor4&   surface_rmatrix,
          Matrix&    surface_emission,
    const Vector&    f_grid,
    const Index&     stokes_dim,
    const Index&     atmosphere_dim,
    const Vector&    rte_los,
    const Numeric&   surface_skin_t,
    const Matrix&    complex_n,
    const Agenda&    blackbody_radiation_agenda,
    const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
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
      surface_specular_R_and_b( ws, surface_rmatrix(0,iv,joker,joker), 
                                surface_emission(iv,joker), Rv, Rh, 
                                f_grid[iv], stokes_dim, surface_skin_t,
                                blackbody_radiation_agenda );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatReflectivity(
          Workspace& ws,
          Matrix&    surface_los,
          Tensor4&   surface_rmatrix,
          Matrix&    surface_emission,
    const Vector&    f_grid,
    const Index&     stokes_dim,
    const Index&     atmosphere_dim,
    const Vector&    rte_los,
    const Numeric&   surface_skin_t,
    const Tensor3&   surface_reflectivity,
    const Agenda&    blackbody_radiation_agenda,
    const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );

  const Index   nf = f_grid.nelem();

  if( surface_reflectivity.nrows() != stokes_dim  &&  
      surface_reflectivity.ncols() != stokes_dim )
    {
      ostringstream os;
      os << "The number of rows and columnss in *surface_reflectivity* must\n"
         << "match *stokes_dim*."
         << "\n stokes_dim : " << stokes_dim 
         << "\n number of rows in *surface_reflectivity* : " 
         << surface_reflectivity.nrows()
         << "\n number of columns in *surface_reflectivity* : " 
         << surface_reflectivity.ncols()
         << "\n";
      throw runtime_error( os.str() );
    }

  if( surface_reflectivity.npages() != nf  &&  
      surface_reflectivity.npages() != 1 )
    {
      ostringstream os;
      os << "The number of pages in *surface_reflectivity* should\n"
         << "match length of *f_grid* or be 1."
         << "\n length of *f_grid* : " << nf 
         << "\n dimension of *surface_reflectivity* : " 
         << surface_reflectivity.npages()
         << "\n";
      throw runtime_error( os.str() );
    }

  out2 << "  Sets variables to model a flat surface\n";

  surface_los.resize( 1, rte_los.nelem() );
  surface_los(0,joker) = rte_los;
  surface_specular_los( surface_los(0, joker), atmosphere_dim );

  surface_emission.resize( nf, stokes_dim );
  surface_rmatrix.resize(1,nf,stokes_dim,stokes_dim);

  Matrix R, IR(stokes_dim,stokes_dim); 

  Vector b;
  blackbody_radiation_agendaExecute( ws, b, surface_skin_t, f_grid, 
                                     blackbody_radiation_agenda ); 
  Vector B(stokes_dim,0);

  for( Index iv=0; iv<nf; iv++ )
    { 
      if( iv == 0  || surface_reflectivity.npages() > 1 )
        { 
          R = surface_reflectivity(iv,joker,joker); 
          for( Index i=0; i<stokes_dim; i++ )
            {
              for( Index j=0; j<stokes_dim; j++ )
                {
                  if( i== j )
                    { IR(i,j) = 1 - R(i,j); }
                  else
                    { IR(i,j) = -R(i,j); }
                }
            }
        }

      surface_rmatrix(0,iv,joker,joker) = R;

      B[0] = b[iv];
      mult( surface_emission(iv,joker), IR, B );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatScalarReflectivity(
          Workspace& ws,
          Matrix&    surface_los,
          Tensor4&   surface_rmatrix,
          Matrix&    surface_emission,
    const Vector&    f_grid,
    const Index&     stokes_dim,
    const Index&     atmosphere_dim,
    const Vector&    rte_los,
    const Numeric&   surface_skin_t,
    const Vector&    surface_scalar_reflectivity,
    const Agenda&    blackbody_radiation_agenda,
    const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 1 );
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
  surface_specular_los( surface_los(0, joker), atmosphere_dim );

  surface_emission.resize( nf, stokes_dim );
  surface_rmatrix.resize(1,nf,stokes_dim,stokes_dim);
  //surface_rmatrix = 0.0;   Not needed when stojkes_dim forced to be 1
  //surface_emission = 0.0;

  Vector b;
  blackbody_radiation_agendaExecute( ws, b, surface_skin_t, f_grid, 
                                     blackbody_radiation_agenda ); 

  Numeric r = 0.0;

  for( Index iv=0; iv<nf; iv++ )
    { 
      if( iv == 0  || surface_scalar_reflectivity.nelem() > 1 )
        { r = surface_scalar_reflectivity[iv]; }

      surface_emission(iv,0) = (1.0-r) * b[iv];
      surface_rmatrix(0,iv,0,0) = r;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceLambertianSimple(
          Workspace& ws,
          Matrix&    surface_los,
          Tensor4&   surface_rmatrix,
          Matrix&    surface_emission,
    const Vector&    f_grid,
    const Index&     stokes_dim,
    const Index&     atmosphere_dim,
    const Vector&    rte_los,
    const Numeric&   surface_skin_t,
    const Vector&    surface_scalar_reflectivity,
    const Agenda&    blackbody_radiation_agenda,
    const Index&     np,
    const Numeric&   za_pos,
    const Verbosity&)
{
  const Index   nf = f_grid.nelem();

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_not_negative( "surface_skin_t", surface_skin_t );
  chk_if_in_range( "za_pos", za_pos, 0, 1 );

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
  const Numeric dza = 90.0 / (Numeric)np;
  const Vector za_lims( 0.0, np+1, dza );

  // surface_los
  for( Index ip=0; ip<np; ip++ )
    { surface_los(ip,0) = za_lims[ip] + za_pos * dza; }

  Vector b;
  blackbody_radiation_agendaExecute( ws, b, surface_skin_t, f_grid, 
                                     blackbody_radiation_agenda ); 

  // Loop frequencies and set remaining values
  //
  Numeric r = 0.0;
  //
  for( Index iv=0; iv<nf; iv++ )
    {
      // Get reflectivity
      if( iv == 0  || surface_scalar_reflectivity.nelem() > 1 )
        { r = surface_scalar_reflectivity[iv]; }

      // surface_rmatrix
      for( Index ip=0; ip<np; ip++ )
        {
          const Numeric w = r * 0.5 * ( cos(2*DEG2RAD*za_lims[ip]) - 
                                         cos(2*DEG2RAD*za_lims[ip+1]) );
          surface_rmatrix(ip,iv,0,0) = w;
        }

      // surface_emission
      surface_emission(iv,0) = (1-r) * b[iv];
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surface_reflectivityFromGriddedField6(
          Tensor3&         surface_reflectivity,
    const Index&           stokes_dim,
    const Vector&          f_grid,
    const Index&           atmosphere_dim,
    const Vector&          lat_grid,
    const Vector&          lat_true,
    const Vector&          lon_true,
    const Vector&          rte_pos,
    const Vector&          rte_los,
    const GriddedField6&   r_field,
    const Verbosity&)
{
  // Basic checks and sizes
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_latlon_true( atmosphere_dim, lat_grid, lat_true, lon_true );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_los( atmosphere_dim, rte_los );
  r_field.checksize_strict();
  //
  const Index nf_in = r_field.data.nvitrines();
  const Index ns2   = r_field.data.nshelves();
  const Index ns1   = r_field.data.nbooks();
  const Index nza   = r_field.data.npages();
  const Index nlat  = r_field.data.nrows();
  const Index nlon  = r_field.data.ncols();
  //
  if( nlat < 2  ||  nlon < 2 )
    {
      ostringstream os;
      os << "The data in *r_field* must span a geographical region. That is,\n"
         << "the latitude and longitude grids must have a length >= 2.";
      throw runtime_error( os.str() );      
    } 
  //
  if( nza < 2 )
    {
      ostringstream os;
      os << "The data in *r_field* must span a range of zenith angles. That\n"
         << "is the zenith angle grid must have a length >= 2.";
      throw runtime_error( os.str() );      

    } 
  if( ns1 < stokes_dim  ||  ns2 < stokes_dim  ||  ns1 > 4  ||  ns2 > 4 )
    {
      ostringstream os;
      os << "The \"Stokes dimensions\" must have a size that is >= "
         << "*stokes_dim* (but not exceeding 4).";
      throw runtime_error( os.str() );      
    } 

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon( lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, 
                                                           lon_true, rte_pos );

  // Interpolate in lat and lon
  Tensor4 r_f_za( nf_in, stokes_dim, stokes_dim, nza );
  {
    GridPos gp_lat, gp_lon;
    gridpos( gp_lat, r_field.get_numeric_grid(4), lat[0] );
    gridpos( gp_lon, r_field.get_numeric_grid(5), lon[0] );
    Vector itw(4);
    interpweights( itw, gp_lat, gp_lon );
    for( Index iv=0; iv<nf_in; iv++ )
      { for( Index iz=0; iz<nza; iz++ )
          { for( Index is1=0; is1<stokes_dim; is1++ )
              { for( Index is2=0; is2<stokes_dim; is2++ )
                  { 
                    r_f_za(iv,is1,is2,iz) = interp( itw, 
                     r_field.data(iv,is1,is2,iz,joker,joker), gp_lat, gp_lon );
  }   }   }   }   }
  
  // Interpolate in incidence angle, cubic if possible
  Tensor3 r_f( nf_in, stokes_dim, stokes_dim );
  Index order = 3;
  if( nza < 4 )
    { order = 1; }
  {
    ArrayOfGridPosPoly   gp(1);
    Matrix               itw(1,order+1);
    Vector               tmp(1);
    gridpos_poly( gp, r_field.get_numeric_grid(3), Vector(1,180-rte_los[0]), 
                                                                       order );
    interpweights( itw, gp );
    //
    for( Index i=0; i<nf_in; i++ )
      { for( Index is1=0; is1<stokes_dim; is1++ )
          { for( Index is2=0; is2<stokes_dim; is2++ )
              { 
                interp( tmp, itw, r_f_za(i,is1,is2,joker), gp );
                r_f(i,is1,is2) = tmp[0];
  }   }   }   }

  // Extract or interpolate in frequency
  //
  if( nf_in == 1 )
    { surface_reflectivity = r_f; }
  else
    {
      const Index nf_out = f_grid.nelem();
      surface_reflectivity.resize( nf_out, stokes_dim, stokes_dim );
      //
      ArrayOfGridPos gp( nf_out );
      Matrix         itw( nf_out, 2 );
      gridpos( gp, r_field.get_numeric_grid(0), f_grid );
      interpweights( itw, gp );
      for( Index is1=0; is1<stokes_dim; is1++ )
        { for( Index is2=0; is2<stokes_dim; is2++ )
            { 
              interp( surface_reflectivity(joker,is1,is2), itw, 
                                       r_f(joker,is1,is2), gp );
    }   }   }     
}



/* Workspace method: Doxygen documentation will be auto-generated */
void surface_scalar_reflectivityFromGriddedField4(
          Vector&          surface_scalar_reflectivity,
    const Index&           stokes_dim,
    const Vector&          f_grid,
    const Index&           atmosphere_dim,
    const Vector&          lat_grid,
    const Vector&          lat_true,
    const Vector&          lon_true,
    const Vector&          rte_pos,
    const Vector&          rte_los,
    const GriddedField4&   r_field,
    const Verbosity&)
{
  // Basic checks and sizes
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 1 );
  chk_latlon_true( atmosphere_dim, lat_grid, lat_true, lon_true );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_los( atmosphere_dim, rte_los );
  r_field.checksize_strict();
  //
  const Index nf_in = r_field.data.nbooks();
  const Index nza   = r_field.data.npages();
  const Index nlat  = r_field.data.nrows();
  const Index nlon  = r_field.data.ncols();
  //
  if( nlat < 2  ||  nlon < 2 )
    {
      ostringstream os;
      os << "The data in *r_field* must span a geographical region. That is,\n"
         << "the latitude and longitude grids must have a length >= 2.";
      throw runtime_error( os.str() );      

    } 
  //
  if( nza < 2 )
    {
      ostringstream os;
      os << "The data in *r_field* must span a range of zenith angles. That\n"
         << "is the zenith angle grid must have a length >= 2.";
      throw runtime_error( os.str() );      
    } 

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon( lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, 
                                                           lon_true, rte_pos );

  // Interpolate in lat and lon
  Matrix r_f_za( nf_in, nza );
  {
    GridPos gp_lat, gp_lon;
    gridpos( gp_lat, r_field.get_numeric_grid(2), lat[0] );
    gridpos( gp_lon, r_field.get_numeric_grid(3), lon[0] );
    Vector itw(4);
    interpweights( itw, gp_lat, gp_lon );
    for( Index iv=0; iv<nf_in; iv++ )
      {
        for( Index iz=0; iz<nza; iz++ )
          { 
            r_f_za(iv,iz) = interp( itw, r_field.data(iv,iz,joker,joker), 
                                                              gp_lat, gp_lon );
          }
      } 

  }    
  
  // Interpolate in incidence angle, cubic if possible
  Vector r_f( nf_in );
  Index order = 3;
  if( nza < 4 )
    { order = 1; }
  {
    ArrayOfGridPosPoly   gp(1);
    Matrix               itw(1,order+1);
    Vector               tmp(1);
    gridpos_poly( gp, r_field.get_numeric_grid(1), Vector(1,180-rte_los[0]), 
                                                                       order );
    interpweights( itw, gp );
    //
    for( Index i=0; i<nf_in; i++ )
      { 
        interp( tmp, itw, r_f_za(i,joker), gp );
        r_f[i] = tmp[0];
      }
  }

  // Extract or interpolate in frequency
  //
  if( nf_in == 1 )
    {
      surface_scalar_reflectivity.resize( 1 );
      surface_scalar_reflectivity[0] = r_f[0];
    }
  else
    {
      const Index nf_out = f_grid.nelem();
      surface_scalar_reflectivity.resize( nf_out );
      //
      ArrayOfGridPos gp( nf_out );
      Matrix         itw( nf_out, 2 );
      gridpos( gp, r_field.get_numeric_grid(0), f_grid );
      interpweights( itw, gp );
      interp( surface_scalar_reflectivity, itw, r_f, gp );
    }     
}


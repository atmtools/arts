/* Copyright (C) 2002-2007 Stefan Buehler <sbuehler@ltu.se>

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

/*!
  \file   special_interp.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-11-14   
  
  \brief  Interpolation routines for special purposes.
  
  This file contains functions connected to interpolation of
  non-general character. The total general interpolation routines are
  found in interpolation.cc. 

  - Interpolation of atmospheric fields and surfaces: 

  These interpolation functions interpolates a field or a surface to a
  set of points, such as the points of a propagation path. That is,
  "blue" interpolation. The functions assume that the grid positions
  are at hand. If several atmospheric fields shall be interpolated,
  the functions to use are *interp_atmfield_gp2itw* and
  *interp_atmfield_by_itw*, where the first function is called once
  and the second is called for each field to be interpolated. If only
  one field shall be interpolated, the function
  *interp_atmfield_by_gp* is a shortcut for calling both functions
  above. The exists an identical set of functions for interpolating
  surfaces, with names where _atmfield_ is replaced with
  _atmsurface_. 

  Possible atmospheric surfaces are *r_geoid*, *z_ground* and one page
  of *z_field*.


  - Conversion of geometric altitudes to pressure values: 

  To convert a geometric altitude to a pressure results in an
  interpolation of the pressure grid, or more exactly the log of the
  *p_grid*. Such functions are placed in this file for that reason.
  These functions have names ending with "2p", for example itw2p.


  - Various functions connected to the atmospheric fields:
  
  See below to find what functions that exist.
*/

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "check_input.h"
#include "special_interp.h"
#include "messages.h"




/*===========================================================================
  === Interpolation functions for atmospheric grids, fields and surfaces
  ===========================================================================*/

//! fix_gridpos_at_boudary
/*!
  This function fixes grid positions on the boundaries of the grids, on which 
  one wants to interpolate. 
  
  The use of the function can be demonstrated using an example, 
  the interpolation of the *pnd_field* on a propagation path inside the 
  cloudbox. The gridpositions calculated in the ppath_agenda are given with 
  respect to the atmospheric grids (p_grid, lat_grid, lon_grid), whereas the 
  *pnd_field* is defined on a subset of these grids, only inside the cloudbox.
  Hence the gridpositions need to be transformed. 
  If a point is exactly at the lower boundary of the cloudbox, 
  the gridposition index can be either cloudbox_limits[0] or
  cloudbox_limits[0]-1, depending on random numerical errors.  To transform 
  the gridpositions, cloudbox_limits[0] is substracted from the original 
  gridposition index, which results in a new gridpostion index of -1. This is 
  of coarse not allowed in the interpolation routines. 
  
  A more precise  example: 
  If a point is exactly on the lower boundary (say p-index = 12) of cloudbox
  the gridpositions can be
              gp.idx = 12, gp.fd[0] = 1, gp.fd[1] = 0
       or     gp.idx = 11, gp.fd[0] = 0, gp.fd[1] = 1
   In the second case gp_cloud.idx will be -1 which gives an error. 

   \param gp An array of gridpositions
   \param grid_size Length of the data grid

   \author Claudia Emde
   \date 2005-05-27
*/
void fix_gridpos_at_boundary(//Input and Output
                             ArrayOfGridPos& gp,
                             const Index grid_size
                             )
{
  // This limit is chosen arbitrarily
  const Numeric epsilon = 1e-6;
  
  for( Index i = 0; i< gp.nelem(); i++)
    {
      if( gp[i].idx == -1 )//&& abs(gp[i].fd[0]-1)<epsilon )
        {
          if (abs(gp[i].fd[0]-1)>epsilon)
            {
              out1 << "  --- WARNING ---, fix_gridpos_at_boundary encountered a strange "
                   << "Gridpos: idx = " << gp[i].idx << ", fd[0] = " << gp[i].fd[0] 
                   << ", fd[1] = " << gp[i].fd[1];
            }
          gp[i].idx += 1;
          gp[i].fd[0] = 0.;
          gp[i].fd[1] = 1.;
        }
      else if (gp[i].idx == grid_size-1)
        {
          if (abs(gp[i].fd[0])>epsilon)
            {
              out1 << "  --- WARNING ---, fix_gridpos_at_boundary encountered a strange "
                   << "Gridpos: idx = " << gp[i].idx << ", fd[0] = " << gp[i].fd[0] 
                   << ", fd[1] = " << gp[i].fd[1];
            }
          gp[i].idx -= 1;
          gp[i].fd[0] = 1.;
          gp[i].fd[1] = 0.;
        }
      if ((gp[i].idx < 0)||(gp[i].idx > grid_size-1))
        {
          ostringstream os;
          os << "Invalid GridPos: idx = " << gp[i].idx << ", fd[0] = " << gp[i].fd[0] << ", fd[1] = " << gp[i].fd[1];
          throw runtime_error(os.str());
        }
      assert(gp[i].idx > -1);
    }
  
}
        

//! interp_atmfield_gp2itw
/*!
    Converts atmospheric grid positions to weights for interpolation of an
    atmospheric field.

    The function is intended for "blue" interpolation, that is, interpolation
    for a set of positions. 

    The output matrix for interpolation weights are resized inside the
    function.

    The input atmospheric grids are checked to be consistent.

    \param   itw                Output: Interpolation weights.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   p_grid             As the WSV with the same name.
    \param   lat_grid           As the WSV with the same name.
    \param   lon_grid           As the WSV with the same name.
    \param   gp_p               Pressure grid positions.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmfield_gp2itw( 
              Matrix&           itw, 
        const Index&            atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  const Index n = gp_p.nelem();

  if( atmosphere_dim == 1 )
    {
      itw.resize(n,2);
      interpweights( itw, gp_p );
    }

  else if( atmosphere_dim == 2 )
    {
      assert( gp_lat.nelem() == n );
      itw.resize(n,4);
      interpweights( itw, gp_p, gp_lat );
    }

  else if( atmosphere_dim == 3 )
    {
      assert( gp_lat.nelem() == n );
      assert( gp_lon.nelem() == n );
      itw.resize(n,8);
      interpweights( itw, gp_p, gp_lat, gp_lon );
    }
}



//! interp_atmfield_by_itw
/*!
    Interpolates an atmospheric field with pre-calculated weights by
    interp_atmfield_gp2itw.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric field is checked to be consistent with the 
    *atmosphere_dim*, *p_grid*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    \param   x                  Output: Values obtained by the interpolation.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   p_grid             As the WSV with the same name.
    \param   lat_grid           As the WSV with the same name.
    \param   lon_grid           As the WSV with the same name.
    \param   x_field            The atmospheric field to be interpolated.
    \param   x_field_name       The name of the field as a string, e.g. 
                                "t_field".
    \param   gp_p               Pressure grid positions.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.
    \param   itw                Interpolation weights from 
                                interp_atmfield_gp2itw.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmfield_by_itw( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        ConstTensor3View        x_field,
        const String&           x_field_name,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw )
{
  // Check consistency between field and grids
  chk_atm_field( x_field_name, x_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  assert( x.nelem() == gp_p.nelem() );

  if( atmosphere_dim == 1 )
    { 
      assert( itw.ncols() == 2 );
      interp( x, itw, x_field(Range(joker),0,0), gp_p ); 
    }

  else if( atmosphere_dim == 2 )
    { 
      assert( itw.ncols() == 4 );
      interp( x, itw, x_field(Range(joker),Range(joker),0), gp_p, gp_lat ); 
    }

  else if( atmosphere_dim == 3 )
    {       
      assert( itw.ncols() == 8 );
      interp( x, itw, x_field, gp_p, gp_lat, gp_lon ); 
    }
}



//! interp_atmfield_by_gp
/*!
    Interpolates an atmospheric field given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric field is checked to be consistent with the 
    *atmosphere_dim*, *p_grid*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    \param   x                  Output: Values obtained by the interpolation.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   p_grid             As the WSV with the same name.
    \param   lat_grid           As the WSV with the same name.
    \param   lon_grid           As the WSV with the same name.
    \param   x_field            The atmospheric field to be interpolated.
    \param   x_field_name       The name of the field as a string, e.g. 
                                "t_field".
    \param   gp_p               Pressure grid positions.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmfield_by_gp( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        ConstTensor3View        x_field,
        const String&           x_field_name,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
  Matrix itw;

  interp_atmfield_gp2itw( itw, atmosphere_dim, p_grid, lat_grid, 
                                              lon_grid, gp_p, gp_lat, gp_lon );

  interp_atmfield_by_itw( x, atmosphere_dim, p_grid, lat_grid, lon_grid,
                            x_field, x_field_name, gp_p, gp_lat, gp_lon, itw );
}



//! interp_atmfield_by_gp
/*!
    As the function above but return-numeric version.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
Numeric interp_atmfield_by_gp( 
        const Index&            atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        ConstTensor3View        x_field,
        const String&           x_field_name,
        const GridPos&          gp_p,
        const GridPos&          gp_lat,
        const GridPos&          gp_lon )
{
  ArrayOfGridPos agp_p(1), agp_lat(0), agp_lon(0);

  gridpos_copy( agp_p[0], gp_p  ); 

  if( atmosphere_dim > 1 )
    {
      agp_lat.resize(1);
      gridpos_copy( agp_lat[0], gp_lat  );
    }

  if( atmosphere_dim > 2 )
    {
      agp_lon.resize(1);
      gridpos_copy( agp_lon[0], gp_lon  );
    }

  Vector x(1);

  interp_atmfield_by_gp( x, atmosphere_dim, p_grid, lat_grid, lon_grid,
                              x_field, x_field_name, agp_p, agp_lat, agp_lon );

  return x[0];
}



//! interp_atmsurface_gp2itw
/*!
    Converts atmospheric grid positions to weights for interpolation of an
    atmospheric surface.

    The function is intended for "blue" interpolation, that is, interpolation
    for a set of positions. 

    The output matrix for interpolation weights are resized inside the
    function.

    The input atmospheric grids are checked to be consistent.

    \param   itw                Output: Interpolation weights.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmsurface_gp2itw( 
              Matrix&           itw, 
        const Index&            atmosphere_dim,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
  if( atmosphere_dim == 1 )
    {
      itw.resize(1,1);
      itw = 1;
    }

  else if( atmosphere_dim == 2 )
    {
      const Index n = gp_lat.nelem();
      itw.resize(n,2);
      interpweights( itw, gp_lat );
    }

  else if( atmosphere_dim == 3 )
    {
      const Index n = gp_lat.nelem();
      assert( n == gp_lon.nelem() );
      itw.resize(n,4);
      interpweights( itw, gp_lat, gp_lon );
    }
}



//! interp_atmsurface_by_itw
/*!
    Interpolates an atmospheric surface with pre-calculated weights by
    interp_atmsurface_gp2itw.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric surface is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    \param   x                  Output: Values obtained by the interpolation.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   lat_grid           As the WSV with the same name.
    \param   lon_grid           As the WSV with the same name.
    \param   x_surface          The atmospheric field to be interpolated.
    \param   x_surface_name     The name of the field as a string, e.g. 
                                "r_geoid".
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.
    \param   itw                Interpolation weights from 
                                interp_atmsurface_gp2itw.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmsurface_by_itw( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        ConstMatrixView         x_surface,
        const String&           x_surface_name,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw )
{
  // Check consistency between field and grids
  chk_atm_surface( x_surface_name, x_surface, atmosphere_dim, lat_grid, 
                                                                    lon_grid );

  if( atmosphere_dim == 1 )
    { 
      assert( itw.ncols() == 1 );
      x = x_surface(0,0); 
    }

  else if( atmosphere_dim == 2 )
    { 
      assert( x.nelem() == gp_lat.nelem() );
      assert( itw.ncols() == 2 );
      interp( x, itw, x_surface(Range(joker),0), gp_lat ); 
    }

  else if( atmosphere_dim == 3 )
    { 
      assert( x.nelem() == gp_lat.nelem() );
      assert( itw.ncols() == 4 );
      interp( x, itw, x_surface, gp_lat, gp_lon ); 
    }
}



//! interp_atmsurface_by_gp
/*!
    Interpolates an atmospheric surface given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric surface is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    \param   x                  Output: Values obtained by the interpolation.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   lat_grid           As the WSV with the same name.
    \param   lon_grid           As the WSV with the same name.
    \param   x_surface          The atmospheric field to be interpolated.
    \param   x_surface_name     The name of the field as a string, e.g. 
                                "t_field".
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmsurface_by_gp( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        ConstMatrixView         x_surface,
        const String&           x_surface_name,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
  Matrix itw;

  interp_atmsurface_gp2itw( itw, atmosphere_dim, gp_lat, gp_lon );

  interp_atmsurface_by_itw( x, atmosphere_dim, lat_grid, lon_grid, 
                              x_surface, x_surface_name, gp_lat, gp_lon, itw );
}



//! interp_atmsurface_by_gp
/*!
    As the function above, but return-numeric version.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
Numeric interp_atmsurface_by_gp( 
        const Index&       atmosphere_dim,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstMatrixView    x_surface,
        const String&      x_surface_name,
        const GridPos&     gp_lat,
        const GridPos&     gp_lon )
{
  ArrayOfGridPos agp_lat(0), agp_lon(0);

  if( atmosphere_dim > 1 )
    {
      agp_lat.resize(1);
      gridpos_copy( agp_lat[0], gp_lat  );
    }

  if( atmosphere_dim > 2 )
    {
      agp_lon.resize(1);
      gridpos_copy( agp_lon[0], gp_lon  );
    }

  Vector x(1);

  interp_atmsurface_by_gp( x, atmosphere_dim, lat_grid, lon_grid, x_surface,
                                            x_surface_name, agp_lat, agp_lon );

  return x[0];
}





/*===========================================================================
  === Conversion altitudes / pressure
  ===========================================================================*/

//! itw2p
/*!
    Converts interpolation weights to pressures.

    The function takes interpolation weights calculated with respect to the 
    vertical dimension, and determines the corresponding pressures. This 
    function can be used when a geometrical altitude is known and the
    pressure for that altitude shall be determined. The interpolation weights
    are then calculated using the geometrical altitudes for the pressure
    surfaces for the position of concern.

    This can be seen as a 1D "blue" interpolation. That means that the number
    of columns of itw shall be 2.

    \param   p_values   Output: Found pressure values.
    \param   p_grid     As the WSV with the same name.
    \param   gp         Altitude grid positions.
    \param   itw        Interpolation weights

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void itw2p(
              VectorView       p_values,
        ConstVectorView        p_grid,
        const ArrayOfGridPos&  gp,
        ConstMatrixView        itw )
{
  assert( itw.ncols() == 2 );
  assert( p_values.nelem() == itw.nrows() );

  // Local variable to store log of the pressure grid:
  Vector logpgrid( p_grid.nelem() );

  transform( logpgrid, log, p_grid );
  
  interp( p_values, itw, logpgrid, gp ); 

  transform( p_values, exp, p_values );
}



//! p2gridpos
/*!
   Calculates grid positions for pressure values.

   This function works as *gridpos*, but is adapted to handle
   pressure grids. The ARTS defintions result in that pressures shall
   not be interpolated directly, it is the log of the pressure that
   shall be interpolated. This means that if some values shall be
   interpolated to some given pressures, the grid positions shall be
   calculated with this function. The interpolation can then be
   performed as usual.

   \param   gp          Output: Grid position Array.
   \param   old_pgrid   The original pressure grid.
   \param   new_pgrid   The new pressure grid.

   \author Patrick Eriksson
   \date   2003-01-20
*/
void p2gridpos(
             ArrayOfGridPos&   gp,
      ConstVectorView          old_pgrid,
      ConstVectorView          new_pgrid )
{
  // Local variable to store log of the pressure grids
  Vector logold( old_pgrid.nelem() );
  Vector lognew( new_pgrid.nelem() );

  transform( logold, log, old_pgrid );
  transform( lognew, log, new_pgrid );
  
  gridpos( gp, logold, lognew );
}

//! p2gridpos_extpol
/*!
   Calculates grid positions for pressure values.

   This function works as *gridpos_extpol*, but is adapted to handle
   pressure grids. The ARTS defintions result in that pressures shall
   not be interpolated directly, it is the log of the pressure that
   shall be interpolated. This means that if some values shall be
   interpolated to some given pressures, the grid positions shall be
   calculated with this function. The interpolation can then be
   performed as usual.

   \param   gp          Output: Grid position Array.
   \param   old_pgrid   The original pressure grid.
   \param   new_pgrid   The new pressure grid.

   \author Patrick Eriksson
   \date   2003-01-20
*/
void p2gridpos_extpol(
             ArrayOfGridPos&   gp,
      ConstVectorView          old_pgrid,
      ConstVectorView          new_pgrid,   
      const Numeric&           extpolfac ) 
{
  // Local variable to store log of the pressure grids
  Vector logold( old_pgrid.nelem() );
  Vector lognew( new_pgrid.nelem() );

  transform( logold, log, old_pgrid );
  transform( lognew, log, new_pgrid );
  
  gridpos_extpol( gp, logold, lognew, extpolfac );
}





/*===========================================================================
  === Various functions connected to the atmospheric fields
  ===========================================================================*/

//! z_at_lat_2d
/*!
    Returns the geomtrical altitudes of *p_grid* for one latitude.

    The latitude is specified by its grid position, in an
    ArrayOfGridPos of length 1. The altitude field (*z_field*) is then
    interpolated to that latitude.

    \param   z          Out: Found altitudes.
    \param   p_grid     As the WSV with the same name.
    \param   lat_grid   As the WSV with the same name.
    \param   z_field    The pressure and latitude part of the WSV with 
                        the same name (that is, the first column).
    \param   gp_lat     Latitude grid position.

    \author Patrick Eriksson 
    \date   2002-11-18
*/
void z_at_lat_2d(
             VectorView   z,
        ConstVectorView   p_grid,
// FIXME only used in assertion
#ifndef NDEBUG
        ConstVectorView   lat_grid,
#else
        ConstVectorView,
#endif
        ConstMatrixView   z_field,
        const GridPos&    gp_lat )
{
  const Index   np = p_grid.nelem();

  assert( z.nelem() == np );
  assert( z_field.nrows() == np );
  assert( z_field.ncols() == lat_grid.nelem() );

  Matrix           z_matrix(np,1);
  ArrayOfGridPos   gp_z(np), agp_lat(1);
  Tensor3          itw(np,1,4);

  gridpos_copy( agp_lat[0], gp_lat );
  gridpos( gp_z, p_grid, p_grid );
  interpweights( itw, gp_z, agp_lat );
  interp( z_matrix, itw, z_field, gp_z, agp_lat );

  z = z_matrix(Range(joker),0);
}



//! z_at_latlon
/*!
    Returns the geomtrical altitudes of *p_grid* for one latitude and
    one longitude.

    The latitude and longitude are specified by their grid position,
    in an ArrayOfGridPos of length 1. The altitude field (*z_field*)
    is then interpolated to that latitude and longitude.

    \param   z          Out: Found altitudes.
    \param   p_grid     As the WSV with the same name.
    \param   lat_grid   As the WSV with the same name.
    \param   lon_grid   As the WSV with the same name.
    \param   z_field    As the WSV with the same name.
    \param   gp_lat     Latitude grid positions.
    \param   gp_lon     Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-12-31
*/
void z_at_latlon(
             VectorView    z,
        ConstVectorView    p_grid,
//FIXME only used in assertion
#ifndef NDEBUG
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
#else
        ConstVectorView,
        ConstVectorView,
#endif
        ConstTensor3View   z_field,
        const GridPos&     gp_lat,
        const GridPos&     gp_lon )
{
  const Index   np = p_grid.nelem();

  assert( z.nelem() == np );
  assert( z_field.npages() == np );
  assert( z_field.nrows() == lat_grid.nelem() );
  assert( z_field.ncols() == lon_grid.nelem() );

  Tensor3          z_tensor(np,1,1);
  ArrayOfGridPos   agp_z(np), agp_lat(1), agp_lon(1);
  Tensor4          itw(np,1,1,8);

  gridpos_copy( agp_lat[0], gp_lat );
  gridpos_copy( agp_lon[0], gp_lon );
  gridpos( agp_z, p_grid, p_grid );
  interpweights( itw, agp_z, agp_lat, agp_lon );

  interp( z_tensor, itw, z_field, agp_z, agp_lat, agp_lon );

  z = z_tensor(Range(joker),0,0);
}

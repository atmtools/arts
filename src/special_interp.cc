/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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
*/

#include <cmath>
#include "check_input.h"
#include "special_interp.h"




//! ArrayOfGridPosPrint
/*!
   Prints a variable of type ArrayOfGridPos to the screen.

   \author Patrick Eriksson
   \date   2002-06-09
*/
void ArrayOfGridPosPrint(
        const ArrayOfGridPos&   x,
        const String&           x_name )
{
  cout << "  *" << x_name <<"*:\n";
  for( Index i=0; i<x.nelem(); i++ )
    cout << "     " << x[i].idx << "  " << x[i].fd[0] << "  " << x[i].fd[1] 
	 << "\n";
}



/*===========================================================================
  === Interpolation functions for atmospheric grids, fields and surfaces
  ===========================================================================*/

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
        const Index&          	atmosphere_dim,
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
        const Index&          	atmosphere_dim,
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
        const Index&          	atmosphere_dim,
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
        const Index&          	atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstTensor3View        x_field,
 	const String&           x_field_name,
        const GridPos&   	gp_p,
        const GridPos&   	gp_lat,
	const GridPos&   	gp_lon )
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
    \param   lat_grid           As the WSV with the same name.
    \param   lon_grid           As the WSV with the same name.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmsurface_gp2itw( 
              Matrix&           itw, 
        const Index&          	atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
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
        const Index&          	atmosphere_dim,
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
    \param   x_field            The atmospheric field to be interpolated.
    \param   x_field_name       The name of the field as a string, e.g. 
                                "t_field".
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmsurface_by_gp( 
              VectorView        x, 
        const Index&          	atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstMatrixView         x_surface,
 	const String&           x_surface_name,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon )
{
  Matrix itw;

  interp_atmsurface_gp2itw( itw, atmosphere_dim, lat_grid, 
 			                            lon_grid, gp_lat, gp_lon );

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



//! itw2p
/*!
    Converts interpolation weights to pressures.

    The function takes interpolation weights calculated with respect to the 
    vertical dimension, and determines the corresponding pressures. This 
    function can be used when a geometrical altitude is known and the
    pressure for that altitude shall be determined. The interpolation weights
    are then calculated using the geometrical altitudes for the pressure
    surfaces for the position of concern.

    This can be seen as a 1D "red" interpolation. That means that the number
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
	const ArrayOfGridPos   gp,
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

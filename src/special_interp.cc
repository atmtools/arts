/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-11-14   
  
  \brief  Interpolation routines for special purposes.
  
  This file contains functions connected to interpolation of
  non-general character. The total general interpolation routines are
  found in interpolation.cc. 

  - Point(s) interpolation of atmospheric/surface fields: 

  These interpolation functions interpolate a atmospheric/surface fields
  to a set of points, such as the points of a propagation path. That is,
  "blue" interpolation. The functions assume that the grid positions
  are at hand. If several atmospheric fields shall be interpolated,
  the functions to use are *interp_atmfield_gp2itw* and
  *interp_atmfield_by_itw*, where the first function is called once
  and the second is called for each field to be interpolated. If only
  one field shall be interpolated, the function
  *interp_atmfield_by_gp* is a shortcut for calling both functions
  above. There exist an identical set of functions for interpolating
  surface-type variables, with names where _atmfield_ is replaced with
  _atmsurface_. 

  Possible surface-type variables are *z_surface* and one page
  of *z_field*.


  - Regridding of atmospheric/surface fields: 

  These functions interpolate from one set of atmospheric grids to a new 
  set of grids.


  - Conversion of geometric altitudes to pressure values: 

  To convert a geometric altitude to a pressure results in an
  interpolation of the pressure grid, or more exactly the log of the
  *p_grid*. Such functions are placed in this file for that reason.
  These functions have names ending with "2p", for example itw2p.


  - Interpolation of Gridded fields of special types:
*/

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "special_interp.h"




/*===========================================================================
  === Point interpolation functions for atmospheric grids and fields
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
    \param   gp_p               Pressure grid positions.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmfield_gp2itw( 
              Matrix&           itw, 
        const Index&            atmosphere_dim,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
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
    \param   x_field            The atmospheric field to be interpolated.
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
        ConstTensor3View        x_field,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw )
{
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
    \param   x_field            The atmospheric field to be interpolated.
    \param   gp_p               Pressure grid positions.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmfield_by_gp( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstTensor3View        x_field,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
  Matrix itw;

  interp_atmfield_gp2itw( itw, atmosphere_dim, gp_p, gp_lat, gp_lon );

  interp_atmfield_by_itw( x, atmosphere_dim, x_field, gp_p, gp_lat, gp_lon, 
                                                                         itw );
}



//! interp_atmfield_by_gp
/*!
    As the function above but return-numeric version.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
Numeric interp_atmfield_by_gp( 
        const Index&            atmosphere_dim,
        ConstTensor3View        x_field,
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

  interp_atmfield_by_gp( x, atmosphere_dim, x_field, agp_p, agp_lat, agp_lon );

  return x[0];
}



//! interp_cloudfield_gp2itw
/*!
    Converts atmospheric a grid position to weights for interpolation of a
    field defined ONLY inside the cloudbox.

    That is, as interp_atmfield_gp2itw, but for cloudbox only variables.

    The input grid position shall be with respect to total grids. If grid
    positions already refer to grid parts inside the cloudbox, you can use 
    interp_atmfield_gp2itw.

    The output grid positions are created by the function, to match the
    cloudbox field, and can be used for later calls of e.g.
    interp_atmfield_by_itw

    \param   itw                Output: Interpolation weights. Vector must be 
                                given correct size before call of function.
    \param   gp_p_out           Output: Pressure cloudbox grid position.
    \param   gp_lat_out         Output: Latitude cloudbox grid positionn.
    \param   gp_lon_out         Output: Longitude cloudbox grid position.
    \param   gp_p_in            Pressure grid position.
    \param   gp_lat_in          Latitude grid position.
    \param   gp_lon_in          Longitude grid position.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   cloudbiox_limits   As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2010-02-12
*/
void interp_cloudfield_gp2itw( 
              VectorView      itw, 
              GridPos&        gp_p_out,
              GridPos&        gp_lat_out,
              GridPos&        gp_lon_out,
        const GridPos&        gp_p_in,
        const GridPos&        gp_lat_in,
        const GridPos&        gp_lon_in,
        const Index&          atmosphere_dim,
        const ArrayOfIndex&   cloudbox_limits )
{
  // Shift grid positions to cloud box grids
  if( atmosphere_dim == 1 )
    {
      gridpos_copy( gp_p_out, gp_p_in ); 
      gp_p_out.idx -= cloudbox_limits[0]; 
      gridpos_upperend_check( gp_p_out, cloudbox_limits[1]-cloudbox_limits[0] );
      assert(itw.nelem() == 2 );
      interpweights( itw, gp_p_out );
    }
  else if( atmosphere_dim == 2 )
    {
      gridpos_copy( gp_p_out, gp_p_in ); 
      gridpos_copy( gp_lat_out, gp_lat_in ); 
      gp_p_out.idx   -= cloudbox_limits[0]; 
      gp_lat_out.idx -= cloudbox_limits[2]; 
      gridpos_upperend_check( gp_p_out,  cloudbox_limits[1]-cloudbox_limits[0]);
      gridpos_upperend_check( gp_lat_out,cloudbox_limits[3]-cloudbox_limits[2]);
      assert(itw.nelem() == 4 );
      interpweights( itw, gp_p_out, gp_lat_out );
    }
  else 
    {
      gridpos_copy( gp_p_out, gp_p_in ); 
      gridpos_copy( gp_lat_out, gp_lat_in ); 
      gridpos_copy( gp_lon_out, gp_lon_in ); 
      gp_p_out.idx   -= cloudbox_limits[0]; 
      gp_lat_out.idx -= cloudbox_limits[2]; 
      gp_lon_out.idx -= cloudbox_limits[4];
      gridpos_upperend_check( gp_p_out,  cloudbox_limits[1]-cloudbox_limits[0]);
      gridpos_upperend_check( gp_lat_out,cloudbox_limits[3]-cloudbox_limits[2]);
      gridpos_upperend_check( gp_lon_out,cloudbox_limits[5]-cloudbox_limits[4]);
      assert(itw.nelem() == 8 );
      interpweights( itw, gp_p_out, gp_lat_out, gp_lon_out );
    }
}



//! interp_atmsurface_gp2itw
/*!
    Converts atmospheric grid positions to weights for interpolation of a
    surface-type variable.

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
    Interpolates a surface-type variable with pre-calculated weights by
    interp_atmsurface_gp2itw.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input surface-type variable is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    \param   x                  Output: Values obtained by the interpolation.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   x_surface          The atmospheric field to be interpolated.
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
        ConstMatrixView         x_surface,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw )
{
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
    Interpolates a surface-type variable given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input surface-type variable is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    \param   x                  Output: Values obtained by the interpolation.
    \param   atmosphere_dim     As the WSV with the same name.
    \param   x_surface          The atmospheric field to be interpolated.
    \param   gp_lat             Latitude grid positions.
    \param   gp_lon             Longitude grid positions.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
void interp_atmsurface_by_gp( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstMatrixView         x_surface,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon )
{
  Matrix itw;

  interp_atmsurface_gp2itw( itw, atmosphere_dim, gp_lat, gp_lon );

  interp_atmsurface_by_itw( x, atmosphere_dim, x_surface, gp_lat, gp_lon, itw );
}



//! interp_atmsurface_by_gp
/*!
    As the function above, but return-numeric version.

    \author Patrick Eriksson 
    \date   2002-11-13
*/
Numeric interp_atmsurface_by_gp( 
        const Index&       atmosphere_dim,
        ConstMatrixView    x_surface,
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

  interp_atmsurface_by_gp( x, atmosphere_dim, x_surface, agp_lat, agp_lon );

  return x[0];
}





/*===========================================================================
  === Regridding
  ===========================================================================*/

//! Regrids an atmospheric field, for precalculated grid positions
/*!
  The function adopts automatically to *atmosphere_dim*. Grid positions not
  used are ignored, i.e. gp_lat is ignored for atmosphere_dim=1 etc.

  \param[out] field_new        Field after interpolation.
  \param[in]  atmosphere_dim   As the WSV with same name.
  \param[in]  field_old        Field to be interpolated.
  \param[in]  gp_p             Pressure grid positions.
  \param[in]  gp_lat           Latitude grid positions.
  \param[in]  gp_lon           Longitude grid positions.

  \author Patrick Eriksson 
  \date   2015-09-09
*/
void regrid_atmfield_by_gp( 
         Tensor3&          field_new, 
   const Index&            atmosphere_dim, 
   ConstTensor3View        field_old, 
   const ArrayOfGridPos&   gp_p,
   const ArrayOfGridPos&   gp_lat,
   const ArrayOfGridPos&   gp_lon )
{
  const Index n1 = gp_p.nelem();

  if( atmosphere_dim == 1 )
    {
      field_new.resize( n1, 1, 1 );
      Matrix itw( n1, 2 );
      interpweights( itw, gp_p );
      interp( field_new(joker,0,0), itw, field_old(joker,0,0), gp_p ); 
    }
  else if( atmosphere_dim == 2 )
    {
      const Index n2 = gp_lat.nelem();
      field_new.resize( n1, n2, 1 );
      Tensor3 itw( n1, n2, 4 );
      interpweights( itw, gp_p, gp_lat );
      interp( field_new(joker,joker,0), itw, field_old(joker,joker,0), gp_p, gp_lat ); 
    }
  else if( atmosphere_dim == 2 )
    {
      const Index n2 = gp_lat.nelem();
      const Index n3 = gp_lon.nelem();
      field_new.resize( n1, n2, n3 );
      Tensor4 itw( n1, n2, n3, 8 );
      interpweights( itw, gp_p, gp_lat, gp_lon );
      interp( field_new, itw, field_old, gp_p, gp_lat, gp_lon ); 
    }
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
    levels for the position of concern.

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

   \param[out]  gp          Output: Grid position Array.
   \param[in]   old_pgrid   The original pressure grid.
   \param[in]   new_pgrid   The new pressure grid.
   \param[in]   extpolfac   Extrapolation factor. Default value is 0.5,
                            which means that extrapolation of half of the
                            last grid distance is allowed.
                            You don't have to specify this.

   \author Patrick Eriksson
   \date   2003-01-20

   \author Stefan Buehler
   \date   2008-03-03
*/
void p2gridpos( ArrayOfGridPos& gp,
                ConstVectorView old_pgrid,
                ConstVectorView new_pgrid,
                const Numeric&  extpolfac )
{
  // Local variable to store log of the pressure grids
  Vector logold( old_pgrid.nelem() );
  Vector lognew( new_pgrid.nelem() );

  transform( logold, log, old_pgrid );
  transform( lognew, log, new_pgrid );

  gridpos( gp, logold, lognew, extpolfac );
}



//! p2gridpos_poly
/*!
 Calculates grid positions for pressure values - higher order interpolation.
 
 This function is similar to p2gridpos, but for higher order interpolation.
 
 \param[out]  gp          Output: Grid position Array.
 \param[in]   old_pgrid   The original pressure grid.
 \param[in]   new_pgrid   The new pressure grid.
 \param[in]   order       Interpolation order (1=linear, 2=quadratic, etc.)
 \param[in]   extpolfac   Extrapolation factor. Default value is 0.5,
                          which means that extrapolation of half of the
                          last grid distance is allowed.
                          You don't have to specify this.
 
 \author Stefan Buehler
 \date   2010-05-03
 */
void p2gridpos_poly( ArrayOfGridPosPoly& gp,
                     ConstVectorView old_pgrid,
                     ConstVectorView new_pgrid,
                     const Index     order,
                     const Numeric&  extpolfac )
{
  // Local variable to store log of the pressure grids
  Vector logold( old_pgrid.nelem() );
  Vector lognew( new_pgrid.nelem() );
  
  transform( logold, log, old_pgrid );
  transform( lognew, log, new_pgrid );
  
  gridpos_poly( gp, logold, lognew, order, extpolfac );
}



//! rte_pos2gridpos
/*!
   Converts a geographical position (rte_pos) to grid positions for p, 
   lat and lon. 

   The function calculates the altitude, latitude and longitude in *rte_pos* to
   matching grid positions. The conversion is straightforwatd for latitude and
   longitude. The altitude shall be converted pressure grid position which
   requires an interpolation of z_field.

   Handles 1D, 2D and 3D (gp_lat and gp_lon untouched if not used).

   Note that the function performs several checks of runtime error type.

   \param   gp_p        Output: Pressure grid position.
   \param   gp_lat      Output: Latitude grid position.
   \param   gp_lon      Output: Longitude grid position.
   \param   atmosphere_dim  As the WSV with the same name.
   \param   p_grid      As the WSV with the same name.
   \param   lat_grid    As the WSV with the same name.
   \param   lon_grid    As the WSV with the same name.
   \param   z_field     As the WSV with the same name.
   \param   rte_pos     As the WSV with the same name.

   \author Patrick Eriksson
   \date   2012-06-22
*/
void rte_pos2gridpos(
         GridPos&     gp_p,
         GridPos&     gp_lat,
         GridPos&     gp_lon,
   const Index&       atmosphere_dim,
   ConstVectorView    p_grid,
   ConstVectorView    lat_grid,
   ConstVectorView    lon_grid,
   ConstTensor3View   z_field,
   ConstVectorView    rte_pos )
{
  chk_rte_pos( atmosphere_dim, rte_pos );

  if( atmosphere_dim == 1 )
    { 
      chk_interpolation_grids( "Altitude interpolation", z_field(joker,0,0), 
                                                                  rte_pos[0] );
      gridpos( gp_p, z_field(joker,0,0), rte_pos[0] ); 
    }
  else
    {
      // Determine z at lat/lon (z_grid) by blue interpolation
      const Index np = p_grid.nelem();
      Vector z_grid( np );
      ArrayOfGridPos agp_z, agp_lat(np);
      //
      gridpos_1to1( agp_z, p_grid );
      //
      chk_interpolation_grids( "Latitude interpolation", lat_grid, rte_pos[1] );
      gridpos( gp_lat, lat_grid, rte_pos[1] );

      if( atmosphere_dim == 2 )
        {
          for( Index i=0; i<np; i++ )
            { agp_lat[i] = gp_lat; }
          Matrix itw( np, 4 );
          interpweights( itw, agp_z, agp_lat );
          interp( z_grid, itw, z_field(joker,joker,0), agp_z, agp_lat );
        }
      else
        {
          chk_interpolation_grids( "Longitude interpolation", lon_grid, 
                                                                  rte_pos[2] );
          gridpos( gp_lon, lon_grid, rte_pos[2] );
          ArrayOfGridPos agp_lon(np);
          for( Index i=0; i<np; i++ )
            { 
              agp_lat[i] = gp_lat;  
              agp_lon[i] = gp_lon; 
            }
          Matrix itw( np, 8 );
          interpweights( itw, agp_z, agp_lat, agp_lon );
          interp( z_grid, itw, z_field, agp_z, agp_lat, agp_lon );

        }

      // And use z_grid to get gp_p (gp_al and gp_lon determined above)
      chk_interpolation_grids( "Altitude interpolation", z_grid, rte_pos[0] );
      gridpos( gp_p, z_grid, rte_pos[0] );
    }
}



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



/*===========================================================================
  === Interpolation of GriddedFields
  ===========================================================================*/

//! complex_n_interp
/*!
    General function for interpolating data of complex n type.

    See documentation of comples_refr_index for format of complex_n-.

    \param   n_real      Out: Real part [nf,nt]
    \param   n_imag      Out: Imaginary part [nf,nt]
    \param   complex_n   Complex refracton index data.
    \param   varname     The name of complex_n to use in error message.
    \param   f_grid      Output frequency grid [nf]
    \param   t_grid      Output temperature grid [nt]

    \author Patrick Eriksson 
    \date   2013-08-16
*/
void complex_n_interp(
         MatrixView       n_real,
         MatrixView       n_imag,
   const GriddedField3&   complex_n,
   const String&          varname,
   ConstVectorView        f_grid,
   ConstVectorView        t_grid )
{
  // Set expected order of grids
  Index gfield_fID    = 0;
  Index gfield_tID    = 1;
  Index gfield_compID = 2;

  // Check of complex_n
  //
  complex_n.checksize_strict();
  //
  chk_griddedfield_gridname( complex_n, gfield_fID,    "Frequency"   );
  chk_griddedfield_gridname( complex_n, gfield_tID,    "Temperature" );
  chk_griddedfield_gridname( complex_n, gfield_compID, "Complex"     );
  //
  if( complex_n.data.ncols() != 2 )
    {
      ostringstream os;
      os << "The data in *" << varname << "* must have exactly two pages. One page "
         << "each\nfor the real and imaginary part of the complex refractive index.";
    } 

  // Frequency and temperature grid sizes
  const Index nf_in = complex_n.data.npages();
  const Index nt_in = complex_n.data.nrows(); 
  const Index nf_out = f_grid.nelem();
  const Index nt_out = t_grid.nelem();

  //Assert size of output variables
  assert( n_real.nrows() == nf_out  &&  n_real.ncols() == nt_out );
  assert( n_imag.nrows() == nf_out  &&  n_imag.ncols() == nt_out );

  const Vector& f_grid_in = complex_n.get_numeric_grid(gfield_fID);
  const Vector& t_grid_in = complex_n.get_numeric_grid(gfield_tID);

  // Expand/interpolate in frequency dimension
  //
  Matrix nrf(nf_out,nt_in), nif(nf_out,nt_in);
  //
  if( nf_in == 1 )
    {
      for( Index i=0; i<nf_out; i++ )
        { 
          nrf(i,joker) = complex_n.data(0,joker,0); 
          nif(i,joker) = complex_n.data(0,joker,1); 
        }
    }
  else
    {
      chk_interpolation_grids( "Frequency interpolation", f_grid_in, f_grid );
      //
      ArrayOfGridPos gp( nf_out );
      Matrix         itw( nf_out, 2 );
      gridpos( gp, f_grid_in, f_grid );
      interpweights( itw, gp );
      for( Index i=0; i<nt_in; i++ )
        { 
          interp( nrf(joker,i), itw, complex_n.data(joker,i,0), gp );
          interp( nif(joker,i), itw, complex_n.data(joker,i,1), gp );
        }
    }

  // Expand/interpolate in temperature dimension
  //
  if( nt_in == 1 )
    {
      for( Index i=0; i<nt_out; i++ )
        { 
          n_real(joker,i) = nrf(joker,0); 
          n_imag(joker,i) = nif(joker,0); 
        }
    }
  else
    {
      chk_interpolation_grids( "Temperature interpolation", t_grid_in, t_grid );
      //
      ArrayOfGridPos gp( nt_out );
      Matrix         itw( nt_out, 2 );
      gridpos( gp, t_grid_in, t_grid );
      interpweights( itw, gp );
      for( Index i=0; i<nf_out; i++ )
        { 
          interp( n_real(i,joker), itw, nrf(i,joker), gp );
          interp( n_imag(i,joker), itw, nif(i,joker), gp );
        }
    }
}

/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
   USA. 
*/

/*!
  \file   cloudbox.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Thu May  23 10:59:55 2002
  
  \brief  Internal functions for scattering calculations.
*/

#include "cloudbox.h"

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <stdexcept>
#include <cmath>

#include "arts.h"
#include "auto_md.h"
#include "make_vector.h"
#include "array.h"
#include "logic.h"
#include "ppath.h"
#include "physics_funcs.h"
#include "array.h"


//! Interpolation of cloud box intensity field
/* 
   See WSM *iyInterpCloudboxField*.

   \param iy                Out: As the WSV with same name.
   \param scat_i_p          In: As the WSV with same name.
   \param scat_i_lat        In: As the WSV with same name.
   \param scat_i_lon        In: As the WSV with same name.
   \param rte_gp_p          In: As the WSV with same name.
   \param rte_gp_lat        In: As the WSV with same name.
   \param rte_gp_lon        In: As the WSV with same name.
   \param rte_los           In: As the WSV with same name.
   \param cloudbox_on       In: As the WSV with same name.
   \param cloudbox_limits   In: As the WSV with same name.
   \param atmosphere_dim    In: As the WSV with same name.
   \param stokes_dim        In: As the WSV with same name.
   \param scat_za_grid      In: As the WSV with same name. 
   \param scat_aa_grid      In: As the WSV with same name. 
   \param f_grid            In: As the WSV with same name.
   \param interpmeth        Interpolation method. Can be "linear" or "cubic".
 
   \author Claudia Emde and Patrick Eriksson
   \date 2004-09-29
*/
void iy_interp_cloudbox_field(
            Matrix&         iy,
      const Tensor7&        scat_i_p,
      const Tensor7&        scat_i_lat,
      const Tensor7&        scat_i_lon,
      const GridPos&        rte_gp_p,
      const GridPos&        rte_gp_lat,
      const GridPos&        rte_gp_lon,
      const Vector&         rte_los,
      const Index&          cloudbox_on,
      const ArrayOfIndex&   cloudbox_limits,
      const Index&          atmosphere_dim,
      const Index&          stokes_dim,
      const Vector&         scat_za_grid,
      const Vector&         scat_aa_grid,
      const Vector&         f_grid,
      const String&         interpmeth )
{
  //--- Check input -----------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  if( !cloudbox_on )
    throw runtime_error( "The cloud box is not activated and no outgoing "
                         "field can be returned." );
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error(
       "*cloudbox_limits* is a vector which contains the upper and lower\n"
       "limit of the cloud for all atmospheric dimensions.\n"
       "So its length must be 2 x *atmosphere_dim*" ); 
  if( scat_za_grid.nelem() == 0 )
    throw runtime_error( "The variable *scat_za_grid* is empty. Are dummy "
                         "values from *cloudboxOff used?" );
  if( !( interpmeth == "linear"  ||  interpmeth == "cubic" ) )
    throw runtime_error( "Unknown interpolation method. Possible choices are "
                                                 "\"linear\" and \"cubic\"." );
  if( interpmeth == "cubic"  &&  atmosphere_dim != 1  )
    throw runtime_error( "Cubic interpolation method is only available for"
                                                     "*atmosphere_dim* = 1." );
  // Size of scat_i_p etc is checked below
  //---------------------------------------------------------------------------


  //--- Determine if at border or inside of cloudbox (or outside!)
  //
  // Let us introduce a number coding for cloud box borders.
  // Borders have the same number as position in *cloudbox_limits*.
  // Innside cloud box is coded as 99, and outside as > 100.
  Index  border  = 999;
  //
  //- Check if at any border
  if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[0] ) )
    { border = 0; }
  else if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[1] ) )
    { border = 1; }
  if( atmosphere_dim == 3  &&  border > 100 )
    {
      if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[2] ) )
        { border = 2; }
      else if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[3] ) )
        { border = 3; }
      else if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[4] ) )
        { border = 4; }
      else if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[5] ) )
        { border = 5; }
    }
  //
  //- Check if inside
  if( border > 100 )
    {
      // Assume inside as it is easiest to detect if outside (can be detected
      // check in one dimension at the time)
      bool inside = true;
      Numeric fgp;

      // Check in pressure dimension
      fgp = fractional_gp( rte_gp_p );
      if( fgp < Numeric(cloudbox_limits[0])  || 
          fgp > Numeric(cloudbox_limits[1]) )
        { inside = false; }

      // Check in lat and lon dimensions
     if( atmosphere_dim == 3  &&  inside )
       {
         fgp = fractional_gp( rte_gp_lat );
         if( fgp < Numeric(cloudbox_limits[2])  || 
             fgp > Numeric(cloudbox_limits[3]) )
           { inside = false; }
         fgp = fractional_gp( rte_gp_lon );
         if( fgp < Numeric(cloudbox_limits[4])  || 
             fgp > Numeric(cloudbox_limits[5]) )
           { inside = false; }
       }

     if( inside )
       { border = 99; }
    }

  // If outside, something is wrong
  if( border > 100 )
    {
      throw runtime_error( 
                 "Given position has been found to be outside the cloudbox." );
    }

  // We are not yet handling points inside the cloud box
  if( border == 99 )
    {
      throw runtime_error( 
               "Observations from inside the cloud box are not yet handled." );
    }

   
  //- Sizes
  const Index   nf  = f_grid.nelem();

  //- Resize *iy*
  iy.resize( nf, stokes_dim );


  // --- 1D ------------------------------------------------------------------
  if( atmosphere_dim == 1 )
    {
      // Use assert to check *scat_i_p* as this variable should to 99% be
      // calculated internally, and thus make it possible to avoid this check.
      assert( is_size( scat_i_p, nf,2,1,1,scat_za_grid.nelem(),1,stokes_dim ));

      out3 << "    Interpolating outgoing field:\n";
      out3 << "       zenith_angle: " << rte_los[0] << "\n";
      if( border )
        out3 << "       top side\n";
      else
        out3 << "       bottom side\n";
      
      // Grid position in *scat_za_grid*
      GridPos gp;
      gridpos( gp, scat_za_grid, rte_los[0] );

      // Corresponding interpolation weights
      Vector itw(2);
      interpweights( itw, gp );

      if( interpmeth == "linear" )
        {
          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                                    scat_i_p( iv, border, 0, 0, joker, 0, is ),
                                      gp );
                }
            }
        }
      else
        {
          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp_cubic( itw, 
                       scat_i_p( iv, border, 0, 0, joker, 0, is ) , rte_los[0],
                                            gp );
                }
            }
        }
    }
  
  // --- 3D ------------------------------------------------------------------
  else
    {
      // Use asserts to check *scat_i_XXX* as these variables should to 99% be
      // calculated internally, and thus make it possible to avoid this check.
      assert ( is_size( scat_i_p, nf, 2, scat_i_p.nshelves(), 
                        scat_i_p.nbooks(), scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(), stokes_dim ));

      assert ( is_size( scat_i_lat, nf, scat_i_lat.nvitrines(), 2, 
                        scat_i_p.nbooks(), scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(), stokes_dim ));
      assert ( is_size( scat_i_lon, nf, scat_i_lat.nvitrines(), 
                        scat_i_p.nshelves(), 2, scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(), stokes_dim ));

      out3 << "    Interpolating outgoing field:\n";
      out3 << "       zenith angle : " << rte_los[0] << "\n";
      out3 << "       azimuth angle: " << rte_los[1]+180 << "\n";

      
      // Scattering angle grid positions
      GridPos gp_za, gp_aa;
      gridpos( gp_za, scat_za_grid, rte_los[0] );
      gridpos( gp_aa, scat_aa_grid, rte_los[1]+180 );

      // Interpolation weights (for 4D "red" interpolation)
      Vector   itw(16);

      // Outgoing from pressure surface
      if( border <= 1 )
        {
          // Lat and lon grid positions with respect to cloud box 
          GridPos cb_gp_lat, cb_gp_lon;
          cb_gp_lat      = rte_gp_lat;
          cb_gp_lon      = rte_gp_lon;
          cb_gp_lat.idx -= cloudbox_limits[2];
          cb_gp_lon.idx -= cloudbox_limits[4]; 
          
          interpweights( itw, cb_gp_lat, cb_gp_lon, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                        scat_i_p( iv, border, joker, joker, joker, joker, is ),
                                      cb_gp_lat, cb_gp_lon, gp_za, gp_aa );
                }
            }
        }

      // Outgoing from latitude surface
      else if( border <= 3 )
        {
          // Pressure and lon grid positions with respect to cloud box 
          GridPos cb_gp_p, cb_gp_lon;
          cb_gp_p        = rte_gp_p;
          cb_gp_lon      = rte_gp_lon;
          cb_gp_p.idx   -= cloudbox_limits[0];
          cb_gp_lon.idx -= cloudbox_limits[4]; 
          
          interpweights( itw, cb_gp_p, cb_gp_lon, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                    scat_i_lat( iv, joker, border-2, joker, joker, joker, is ),
                                      cb_gp_p, cb_gp_lon, gp_za, gp_aa );
                }
            }
        }

      // Outgoing from longitude surface
      else
        {
          // Pressure and lat grid positions with respect to cloud box 
          GridPos cb_gp_p, cb_gp_lat;
          cb_gp_p        = rte_gp_p;
          cb_gp_lat      = rte_gp_lat;
          cb_gp_p.idx   -= cloudbox_limits[0]; 
          cb_gp_lat.idx -= cloudbox_limits[2];
          
          interpweights( itw, cb_gp_p, cb_gp_lat, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                    scat_i_lon( iv, joker, joker, border-4, joker, joker, is ),
                                      cb_gp_p, cb_gp_lat, gp_za, gp_aa );
                }
            }
        }
    }
}


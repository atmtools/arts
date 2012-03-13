/* Copyright (C) 2002-2008 Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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
  ===  File description 
  ===========================================================================*/

/*!
  \file   m_ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions releated to propagation paths variables.

  The file includes special functions to set the sensor position and LOS,
  and functions for calculation of propagation paths.

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
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "xml_io.h"
#include "refraction.h"
#include "m_general.h"

extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;
extern const Numeric EARTH_GRAV_CONST;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCalc(Workspace&            ws,
               // WS Output:
               Ppath&                ppath,
               // WS Input:
               const Agenda&         ppath_step_agenda,
               const Index&          basics_checked,
               const Index&          atmosphere_dim,
               const Vector&         p_grid,
               const Vector&         lat_grid,
               const Vector&         lon_grid,
               const Tensor3&        z_field,
               const Vector&         refellipsoid,
               const Matrix&         z_surface,
               const Index&          cloudbox_on, 
               const Index&          cloudbox_checked,
               const ArrayOfIndex&   cloudbox_limits,
               const Vector&         rte_pos,
               const Vector&         rte_los,
               const Verbosity&      verbosity)
{
  // Basics and cloudbox OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, 
              p_grid, lat_grid, lon_grid, z_field, refellipsoid, z_surface, 
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1,
              verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepGeometric(// WS Output:
                         Ppath&           ppath_step,
                         // WS Input:
                         const Index&     atmosphere_dim,
                         const Vector&    lat_grid,
                         const Vector&    lon_grid,
                         const Tensor3&   z_field,
                         const Vector&    refellipsoid,
                         const Matrix&    z_surface,
                         const Numeric&   ppath_lmax,
                         const Verbosity&)
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here.

  if( atmosphere_dim == 1 )
    { ppath_step_geom_1d( ppath_step, z_field(joker,0,0), 
                          refellipsoid, z_surface(0,0), ppath_lmax ); }

  else if( atmosphere_dim == 2 )
    { ppath_step_geom_2d( ppath_step, lat_grid,
                          z_field(joker,joker,0), refellipsoid, 
                          z_surface(joker,0), ppath_lmax ); }


  else if( atmosphere_dim == 3 )
    { ppath_step_geom_3d( ppath_step, lat_grid, lon_grid,
                          z_field, refellipsoid, z_surface, ppath_lmax ); }

  else
    { throw runtime_error( "The atmospheric dimensionality must be 1-3." ); }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepGeometricUpDown(// WS Output:
                         Ppath&           ppath_step,
                         // WS Input:
                         const Index&     atmosphere_dim,
                         const Tensor3&   z_field,
                         const Vector&    refellipsoid,
                         const Matrix&    z_surface,
                         const Index&     cloudbox_on, 
                    const ArrayOfIndex&   cloudbox_limits,
                         const Numeric&   ppath_lmax,
                         const Verbosity&)
{
  // Check if OK to use function
  if( atmosphere_dim > 1 )
    {
      throw runtime_error( "This method can only be used for 1D atmospheres." );
    }
  if( ppath_lmax >= 0 )
    {
      throw runtime_error( "This method can only be used for ppath_lmax < 0." );
    }
  // Zenith angle checked in sub-function

  ppath_geom_updown_1d( ppath_step, z_field(joker,0,0), refellipsoid, 
                        z_surface(0,0), cloudbox_on, cloudbox_limits );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepRefractionEuler(Workspace&  ws,
                               // WS Output:
                                     Ppath&      ppath_step,
                               // WS Input:
                               const Agenda&     refr_index_agenda,
                               const Index&      atmosphere_dim,
                               const Vector&     p_grid,
                               const Vector&     lat_grid,
                               const Vector&     lon_grid,
                               const Tensor3&    z_field,
                               const Tensor3&    t_field,
                               const Tensor4&    vmr_field,
                               const Vector&     refellipsoid,
                               const Matrix&     z_surface,
                               const Numeric&    ppath_lmax,
                               const Numeric&    ppath_lraytrace,
                               const Verbosity&)
{
  // Input checks here would be rather costly as this function is called
  // many times. 

  assert( ppath_lraytrace > 0 );

  if( atmosphere_dim == 1 )
    { 
      ppath_step_refr_1d( ws, ppath_step, p_grid, z_field(joker,0,0), 
                          t_field(joker,0,0), vmr_field(joker,joker,0,0), 
                          refellipsoid, z_surface(0,0), ppath_lmax, 
                          refr_index_agenda, "linear_euler", ppath_lraytrace );
    }
  else if( atmosphere_dim == 2 )
    { 
      ppath_step_refr_2d( ws, ppath_step, p_grid, lat_grid, 
                          z_field(joker,joker,0), t_field(joker,joker,0), 
                          vmr_field(joker, joker,joker,0), refellipsoid, 
                          z_surface(joker,0), ppath_lmax, refr_index_agenda,
                          "linear_euler", ppath_lraytrace ); 
    }
  else if( atmosphere_dim == 3 )
    { 
      ppath_step_refr_3d( ws, ppath_step, p_grid, lat_grid, lon_grid, z_field, 
                          t_field, vmr_field, refellipsoid, z_surface, 
                          ppath_lmax, refr_index_agenda,
                          "linear_euler", ppath_lraytrace ); 
    }
  else
    { throw runtime_error( "The atmospheric dimensionality must be 1-3." ); }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losSet(// WS Output:
                Vector&          rte_los,
                // WS Input:
                const Index&     atmosphere_dim,
                // Control Parameters:
                const Numeric&   za,
                const Numeric&   aa,
                const Verbosity&)
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  if( atmosphere_dim == 1 )
    { rte_los.resize(1); }
  else
    {
      rte_los.resize(2);
      rte_los[1] = aa;
    }
  rte_los[0] = za;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void rte_posSet(// WS Output:
                Vector&          rte_pos,
                // WS Input:
                const Index&     atmosphere_dim,
                // Control Parameters:
                const Numeric&   z,
                const Numeric&   lat,
                const Numeric&   lon,
                const Verbosity&)
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  rte_pos.resize(atmosphere_dim);
  rte_pos[0] = z;
  if( atmosphere_dim >= 2 )
    { rte_pos[1] = lat; }
  if( atmosphere_dim == 3 )
    { rte_pos[2] = lon; }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZaRefr1D(Workspace&          ws,
                          // WS Generic Output:
                          Vector&             za_vector,
                          // WS Input:
                          const Agenda&       refr_index_agenda,
                          const Matrix&       sensor_pos,
                          const Vector&       p_grid,
                          const Tensor3&      t_field,
                          const Tensor3&      z_field,
                          const Tensor4&      vmr_field,
                          const Vector&       refellipsoid,
                          const Index&        atmosphere_dim,
                          // WS Generic Input:
                          const Vector&       ztan_vector,
                          const Verbosity&)
{
  if( atmosphere_dim != 1 ) {
    throw runtime_error( "The function can only be used for 1D atmospheres." );
  }

  if( ztan_vector.nelem() != sensor_pos.nrows() ) {
    ostringstream os;
    os << "The number of altitudes in true tangent altitude vector must\n"
       << "match the number of positions in *sensor_pos*.";
    throw runtime_error( os.str() );
  }

  // Set za_vector's size equal to ztan_vector
  za_vector.resize( ztan_vector.nelem() );

  // Define refraction variables
  Numeric   refr_index;

  // Calculate refractive index for the tangential altitudes
  for( Index i=0; i<ztan_vector.nelem(); i++ ) {
    get_refr_index_1d( ws, refr_index, refr_index_agenda, p_grid, 
                       refellipsoid[0], z_field(joker,0,0), t_field(joker,0,0), 
                       vmr_field(joker,joker,0,0), 
                       ztan_vector[i] + refellipsoid[0] );

    // Calculate zenith angle
    za_vector[i] = 180 - RAD2DEG* asin( refr_index * 
                                        (refellipsoid[0] + ztan_vector[i]) / 
                                        (refellipsoid[0] + sensor_pos(i,0)) );
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZa1D(// WS Generic Output:
                      Vector&             za_vector,
                      // WS Input:
                      const Matrix&       sensor_pos,
                      const Vector&       refellipsoid,
                      const Index&        atmosphere_dim,
                      // WS Generic Input:
                      const Vector&       ztan_vector,
                      const Verbosity&)
{
  if( atmosphere_dim != 1 ) {
    throw runtime_error( "The function can only be used for 1D atmospheres." );
  }

  const Index   npos = sensor_pos.nrows();

  if( ztan_vector.nelem() != npos )
    {
      ostringstream os;
      os << "The number of altitudes in the geometric tangent altitude vector\n"
         << "must match the number of positions in *sensor_pos*.";
      throw runtime_error( os.str() );
    }

  za_vector.resize( npos );

  for( Index i=0; i<npos; i++ )
    { za_vector[i] = geompath_za_at_r( refellipsoid[0] + ztan_vector[i], 100,
                                       refellipsoid[0] + sensor_pos(i,0) ); }
}



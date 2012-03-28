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
#include "geodetic.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "xml_io.h"
#include "refraction.h"
#include "m_general.h"

extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCalc(      
         Workspace&      ws,
         Ppath&          ppath,
         Vector&         rte_los,
   const Agenda&         ppath_agenda,
   const Index&          basics_checked,
   const Tensor3&        t_field,
   const Tensor3&        z_field,
   const Tensor4&        vmr_field,
   const Index&          cloudbox_on, 
   const Index&          cloudbox_checked,
   const Index&          ppath_inside_cloudbox_do,
   const Index&          mblock_index,
   const Vector&         rte_pos,
   const Verbosity&  )
{
  //--- Check input -----------------------------------------------------------
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  ppath_agendaExecute( ws, ppath, rte_los, rte_pos, cloudbox_on, 
                       ppath_inside_cloudbox_do, mblock_index, 
                       t_field, z_field, vmr_field, ppath_agenda );
}


void line_circle_intersect(
         Numeric&   x,
         Numeric&   z,
   const Numeric&   xl,
   const Numeric&   zl,
   const Numeric&   dx,
   const Numeric&   dz,
   const Numeric&   xc,
   const Numeric&   zc,
   const Numeric&   r )
{
  const Numeric a = dx*dx + dz*dz;
  const Numeric b = 2*( dx*(xl-xc) + dz*(zl-zc) );
  const Numeric c = xc*xc + zc*zc + xl*xl + zl*zl - 2*(xc*xl + zc*zl) - r*r;
  
  Numeric d = b*b - 4*a*c;
  assert( d > 0 );
  
  const Numeric a2 = 2*a;
  const Numeric b2 = -b / a2;
  const Numeric e  = sqrt( d ) / a2;

  const Numeric l1 = b2 + e;
  const Numeric l2 = b2 + e;

  Numeric l;
  if( l1 < 0 )
    { l = l2; }
  else  if( l2 < 0 )
    { l = l1; }
  else
    { l = min(l1,l2); assert( l>=0 ); }

  x = xl + l*dx;
  z = zl + l*dz;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathFromRtePos2(      
         Workspace&      ws,
         Ppath&          ppath,
         Vector&         rte_los,
   const Agenda&         ppath_step_agenda,
   const Index&          basics_checked,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        t_field,
   const Tensor3&        z_field,
   const Tensor4&        vmr_field,
   const Vector&         refellipsoid,
   const Matrix&         z_surface,
   const Vector&         rte_pos,
   const Vector&         rte_pos2,
   const Verbosity&      verbosity )
{
  //--- Check input -----------------------------------------------------------
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  chk_rte_pos( atmosphere_dim, rte_pos, 0 );
  chk_rte_pos( atmosphere_dim, rte_pos2, 1 );

  // Radius of rte_pos and rte_pos2
  const Numeric r1  = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                       lon_grid, rte_pos ) + rte_pos[0];
  const Numeric r2 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                       lon_grid, rte_pos2 ) + rte_pos2[0];

  // LOS from rte_pos to rte_pos2
  Vector rte_los_geom;
  rte_losGeometricFromRtePosToRtePos2( rte_los_geom, atmosphere_dim, lat_grid, 
                       lon_grid, refellipsoid, rte_pos, rte_pos2, verbosity );

  // Geometric distance between rte_pos and rte_pos2
  Numeric l12, lat1=0;
  if( atmosphere_dim <= 2 )
    { 
      if( atmosphere_dim == 2 )
        { lat1 = rte_pos[1]; }
      distance2D( l12, r1, lat1, r2, rte_pos2[1] ); 
    } 
  else
    { distance3D( l12, r1,  rte_pos[1],  rte_pos[2], 
                       r2, rte_pos2[1], rte_pos2[2] ); }
  
  // Iterate until convergence for rte_los
  //
  rte_los     = rte_los_geom;
  Index ready = false, counter=0;
  //
  while( !ready )
    {
      // Path for present rte_los (no cloudbox!)
      ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
                  lat_grid, lon_grid, t_field, z_field, vmr_field, 
                  refellipsoid, z_surface, 0, ArrayOfIndex(0), 
                  rte_pos, rte_los, 0, verbosity );

      // Find the point closest to rte_pos2, on the side towards rte_pos
      // We do this by looking at the distance to rte_pos (should be as close
      // to l12 as possible, but not exceed it)
      Numeric lip = 99e99;
      Index ip = ppath.np;
      //
      while( lip >= l12  &&  ip > 0 )
        {
          ip--;
          
          // Distance between rte_pos and ppath point
          if( atmosphere_dim <= 2 )
            { distance2D( lip, r1, lat1, ppath.r[ip], ppath.pos(ip,1) ); }
          else
            { distance3D( lip, r1, rte_pos[1], rte_pos[2], 
                          ppath.r[ip], ppath.pos(ip,1), ppath.pos(ip,2) ); }
        }

      // Estimate ppath at the distance of l12, and the distance for that point
      // to rte_pos2 (dd)
      Vector  posc;
      Numeric dd;
      {
        if( atmosphere_dim <= 2 )
          { 
            Numeric xip, zip, dxip, dzip, x1, z1, xc, zc, rc, latc;
            poslos2cart( xip, zip, dxip, dzip, ppath.r[ip], 
                         ppath.pos(ip,1), ppath.los(ip,0) );
            pol2cart( x1, z1, r1, lat1 );
            line_circle_intersect( xc, zc, xip, zip, dxip, dzip, x1, z1, l12 );
            cart2pol( rc, latc, xc, zc, ppath.pos(ip,1), ppath.los(ip,0) );
            posc.resize(2);
            posc[1] = latc;
            posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                           lat_grid, lon_grid, posc );
            distance2D( dd, rc, latc, r2, rte_pos2[1] );
            
            /*
            const Numeric dl = l12 - lip;
            cout << "dl = " << dl << endl;
            Numeric rl, latl, xip, zip, dxip, dzip;
            poslos2cart( xip, zip, dxip, dzip, ppath.r[ip], 
                         ppath.pos(ip,1), ppath.los(ip,0) );
            cart2pol( rl, latl, xip+dl*dxip, zip+dl*dzip,
                      ppath.pos(ip,1), ppath.los(ip,0) );
            posc.resize(2);
            posc[1] = latl;
            posc[0] = rl - pos2refell_r( atmosphere_dim, refellipsoid, 
                                           lat_grid, lon_grid, posc );
            cout << "zl = " << posc[0]/1e3 << endl;
            cout << "latl= " << latl << endl;
            distance2D( dd, rl, latl, r2, rte_pos2[1] );
            */
          }
        else
          { 
            const Numeric dl = l12 - lip;
            Numeric rl, latl, lonl, xip, yip, zip, dxip, dyip, dzip;
            poslos2cart( xip, yip, zip, dxip, dyip, dzip, ppath.r[ip], 
                         ppath.pos(ip,1), ppath.pos(ip,2), 
                         ppath.los(ip,0), ppath.los(ip,1) );
            cart2sph( rl, latl, lonl, xip+dl*dxip, yip+dl*dyip, zip+dl*dzip,
                      ppath.pos(ip,1), ppath.pos(ip,2), 
                      ppath.los(ip,0), ppath.los(ip,1) );
            posc.resize(3);
            posc[1] = latl;
            posc[2] = lonl;
            posc[0] = rl - pos2refell_r( atmosphere_dim, refellipsoid, 
                                           lat_grid, lon_grid, posc );
            distance3D( dd, rl, latl, lonl, r2, rte_pos2[1], rte_pos2[2] );
          }
      }

      // Converged?
      if( dd < 1e-3 )
        { ready = true; }

      // Otherwise, calculate correction for rte_los
      else
        {
          Vector los;
          rte_losGeometricFromRtePosToRtePos2( los, atmosphere_dim, lat_grid, 
                        lon_grid, refellipsoid, rte_pos, posc, verbosity );
          rte_los[0] -= los[0] - rte_los_geom[0];
          if( atmosphere_dim == 3 )
            { rte_los[1] -= los[1] - rte_los_geom[1]; }

          cout << rte_los[0]-rte_los_geom[0] << " " << los[0]-rte_los_geom[0]
               << endl;
        }

      // Brake if not finding a solution
      counter++;
      if( counter >= 30 )
        throw runtime_error( "Sorry, but the algorithm is not converging. "
                             "Don't know what to do!" );
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void ppathStepByStep(
         Workspace&      ws,
         Ppath&          ppath,
   const Agenda&         ppath_step_agenda,
   const Index&          ppath_inside_cloudbox_do,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        t_field,
   const Tensor3&        z_field,
   const Tensor4&        vmr_field,
   const Vector&         refellipsoid,
   const Matrix&         z_surface,
   const Index&          cloudbox_on, 
   const ArrayOfIndex&   cloudbox_limits,
   const Vector&         rte_pos,
   const Vector&         rte_los,
   const Verbosity&      verbosity)
{
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid, 
              lon_grid, t_field, z_field, vmr_field, refellipsoid, z_surface, 
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 
              ppath_inside_cloudbox_do, verbosity );
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
  
  // A call with background set, just wants to obtain the refractive index for
  // complete ppaths consistent of a single point.
  if( !ppath_what_background( ppath_step ) )
    {
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
        { throw runtime_error( "The atmospheric dimensionality must be 1-3." );}
    }

  else
    { 
      assert( ppath_step.np == 1 );
      ppath_step.nreal[0] = 1;
    }
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

  // A call with background set, just wants to obtain the refractive index for
  // complete ppaths consistent of a single point.
  if( !ppath_what_background( ppath_step ) )
    {
      if( atmosphere_dim == 1 )
        { 
          ppath_step_refr_1d( ws, ppath_step, p_grid, z_field(joker,0,0), 
                              t_field(joker,0,0), vmr_field(joker,joker,0,0), 
                              refellipsoid, z_surface(0,0), ppath_lmax, 
                              refr_index_agenda, "linear_euler", 
                              ppath_lraytrace );
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
          ppath_step_refr_3d( ws, ppath_step, p_grid, lat_grid, lon_grid, 
                              z_field, t_field, vmr_field, refellipsoid, 
                              z_surface, ppath_lmax, refr_index_agenda,
                              "linear_euler", ppath_lraytrace ); 
        }
      else
        { throw runtime_error( "The atmospheric dimensionality must be 1-3." );}
    }

  else
    { 
      assert( ppath_step.np == 1 );
      if( atmosphere_dim == 1 )
        { get_refr_index_1d( ws, ppath_step.nreal[0], refr_index_agenda, 
                             p_grid, refellipsoid, z_field(joker,0,0), 
                             t_field(joker,0,0), vmr_field(joker,joker,0,0), 
                             ppath_step.r[0] ); 
        }
      else if( atmosphere_dim == 2 )
        { get_refr_index_2d( ws, ppath_step.nreal[0], refr_index_agenda, 
                             p_grid, lat_grid, refellipsoid, 
                             z_field(joker,joker,0), t_field(joker,joker,0), 
                             vmr_field(joker, joker,joker,0),
                             ppath_step.r[0], ppath_step.pos(0,1) ); 
        }
      else
        { get_refr_index_3d( ws, ppath_step.nreal[0], refr_index_agenda, 
                             p_grid, lat_grid, lon_grid, refellipsoid, 
                             z_field, t_field, vmr_field, ppath_step.r[0], 
                             ppath_step.pos(0,1), ppath_step.pos(0,2) ); 
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void PrintTangentPoint(
   const Ppath&     ppath,
   const Index&     level,
   const Verbosity& verbosity)
{
  // Find lowest z and its index
  Numeric zmin = 99e99;
  Index   imin = -1;
  for( Index i=0; i<ppath.np; i++ )
    {
      if( ppath.pos(i,0) < zmin )
        {
          zmin = ppath.pos(i,0);
          imin = i;
        }
    }

  ostringstream os;

  if( imin == 0  ||  imin == ppath.np-1 )
    {
      os << "Lowest altitude found at the end of the propagation path.\n"
         << "This indicates that the tangent point is either above the\n"
         << "top-of-the-atmosphere or below the planet's surface.";
    }

  else
    {
      os << "   z: " << ppath.pos(imin,0)/1e3 << " km\n" 
         << " lat: " << ppath.pos(imin,1) << " deg";
        if( ppath.pos.ncols() == 3 )
          os << "\n lon: " << ppath.pos(imin,2) << " deg";
    }

  CREATE_OUTS
  SWITCH_OUTPUT (level, os.str ());  
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
void rte_losGeometricFromRtePosToRtePos2(
         Vector&         rte_los,
   const Index&          atmosphere_dim,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Vector&         refellipsoid,
   const Vector&         rte_pos,
   const Vector&         rte_pos2,
   const Verbosity& )
{
  // Check input
  chk_rte_pos( atmosphere_dim, rte_pos, 0 );
  chk_rte_pos( atmosphere_dim, rte_pos2, 1 );

  // Radius of rte_pos and rte_pos2
  const Numeric r1  = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                       lon_grid, rte_pos ) + rte_pos[0];
  const Numeric r2 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                       lon_grid, rte_pos2 ) + rte_pos2[0];

  // Remaining polar and cartesian coordinates of rte_pos
  Numeric lat1, lon1=0, x1, y1=0, z1;
  // Cartesian coordinates of rte_pos2
  Numeric x2, y2=0, z2;
  //
  if( atmosphere_dim == 1 )
    {
      // Latitude distance implicitly checked by chk_rte_pos 
      lat1 = 0;
      pol2cart( x1, z1, r1, lat1 );
      pol2cart( x2, z2, r2, rte_pos2[1] );
    }
  else if( atmosphere_dim == 2 )
    {
      lat1 = rte_pos[1];
      pol2cart( x1, z1, r1, lat1 );
      pol2cart( x2, z2, r2, rte_pos2[1] );
    }
  else 
    {
      lat1 = rte_pos[1];
      lon1 = rte_pos[2];
      sph2cart( x1, y1, z1, r1, lat1, lon1 );
      sph2cart( x2, y2, z2, r2, rte_pos2[1], rte_pos2[2] );
    }

  // Geometrical LOS to transmitter
  Numeric za, aa;
  //
  los2xyz( za, aa, r1, lat1, lon1, x1, y1, z1, x2, y2, z2 );
  //
  if( atmosphere_dim == 3 )
    { 
      rte_los.resize(2); 
      rte_los[0] = za; 
      rte_los[1] = aa; 
    }
  else 
    { 
      rte_los.resize(1); 
      rte_los[0] = za; 
      if( atmosphere_dim == 2  && aa < 0 ) // Should 2D-za be negative?
        { rte_los[0] = -za; }
    }
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
  for( Index i=0; i<ztan_vector.nelem(); i++ ) 
    {
      get_refr_index_1d( ws, refr_index, refr_index_agenda, p_grid, 
                         refellipsoid[0], z_field(joker,0,0), 
                         t_field(joker,0,0), vmr_field(joker,joker,0,0), 
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



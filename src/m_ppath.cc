/* Copyright (C) 2002-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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
#include "lin_alg.h"
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
    const Agenda&         ppath_agenda,
    const Numeric&        ppath_lraytrace,
    const Index&          basics_checked,
    const Tensor3&        t_field,
    const Tensor3&        z_field,
    const Tensor4&        vmr_field,
    const Tensor3&        edensity_field,
    const Vector&         f_grid,
    const Index&          cloudbox_on, 
    const Index&          cloudbox_checked,
    const Index&          ppath_inside_cloudbox_do,
    const Vector&         rte_pos,
    const Vector&         rte_los,
    const Vector&         rte_pos2,
    const Verbosity&  )
{
  //--- Check input -----------------------------------------------------------
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  ppath_agendaExecute( ws, ppath, ppath_lraytrace, rte_pos, rte_los, rte_pos2,
                       cloudbox_on, ppath_inside_cloudbox_do, t_field, z_field,
                       vmr_field, edensity_field, f_grid, ppath_agenda );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ppathFromRtePos2V2(      
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
    const Tensor3&        edensity_field,
    const Vector&         f_grid,
    const Vector&         refellipsoid,
    const Matrix&         z_surface,
    const Vector&         rte_pos,
    const Vector&         rte_pos2,
    const Numeric&        ppath_lraytrace,
    const Verbosity&      verbosity )
{
  //--- Check input -----------------------------------------------------------
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_pos( atmosphere_dim, rte_pos2, true );
  //---------------------------------------------------------------------------

 
  // Geometric LOS from rte_pos to rte_pos2
  Vector rte_los_geom;
  rte_losGeometricFromRtePosToRtePos2( rte_los_geom, atmosphere_dim, lat_grid, 
                        lon_grid, refellipsoid, rte_pos, rte_pos2, verbosity );

  // Radius of rte_pos and rte_pos2
  const Numeric r1 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                              lon_grid, rte_pos ) + rte_pos[0];
  const Numeric r2 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                            lon_grid, rte_pos2 ) + rte_pos2[0];

  // Geometric distance between rte_pos and rte_pos2, effective 2D-lat for
  // rte_pos and and Cartesian coordinates of rte_pos:
  Numeric l12, lat1=0, x1, y1=0, z1;
  if( atmosphere_dim <= 2 )
    { 
      if( atmosphere_dim == 2 )
        { lat1 = rte_pos[1]; }
      distance2D( l12, r1, lat1, r2, rte_pos2[1] ); 
      pol2cart( x1, z1, r1, lat1 );
    } 
  else
    { 
      distance3D( l12, r1, rte_pos[1],  rte_pos[2], 
                       r2, rte_pos2[1], rte_pos2[2] ); 
      sph2cart( x1, y1, z1, r1, rte_pos[1], rte_pos[2] );
    }

  // Add handling of 2D with negative za !!!!!!

  // Define remaining variables used in the while-loop below
  //
  // Basic bookkeeping variables
  Numeric za_upp_limit = 180;
  Numeric za_low_limit = 0;
  //
  // Various variables associated with the ppath, and the point of the path
  // closest to the transmitter
  Ppath   ppt;                // "Test ppath"
  Index   ip  = -999;         // Index of closest ppath point
  Numeric dxip, dyip=0, dzip; // Cartesian LOS of the closest ppath point
  //
  // Data for the intersection of the l12-sphere
  Vector  posc(max( Index(2),atmosphere_dim));
  Numeric rc, xc, yc=0, zc;   

  CREATE_OUT2;
  CREATE_OUT3;

  Vector t_za(99), t_dza(99);
  Index  it = -1;

  const Numeric za_accuracy = 1e-5;

  // Keep trying until ready or ground intersetion determined
  //
  bool  ready  = false;
  bool  ground = false;
  bool  failed = false;
  Index ntries = 0;
  //
  while( !ready && !ground )
    {
      out2 << "    Trying rte_los: [" << rte_los << "]\n";

      // Path for present rte_los (no cloudbox!)
      ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
                  lon_grid, t_field, z_field, vmr_field, edensity_field,
                  f_grid, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
                  rte_pos, rte_los, ppath_lraytrace, 0, verbosity );

      // Find the point closest to rte_pos2, on the side towards rte_pos. 
      // We do this by looking at the distance to rte_pos, that should be 
      // as close to l12 as possible, but not exceed it.
      Numeric lip = 99e99;
      ip = ppt.np;
      //
      while( lip >= l12  &&  ip > 0 )
        {
          ip--;
          if( atmosphere_dim <= 2 )
            { distance2D( lip, r1, lat1, ppt.r[ip], ppt.pos(ip,1) ); }
          else
            { distance3D( lip, r1, rte_pos[1], rte_pos[2], 
                          ppt.r[ip], ppt.pos(ip,1), ppt.pos(ip,2) ); }
        }

      Numeric za_new, daa=0;

      // Surface intersection:
      // Not OK if the ground position is too far from rte_pos2. 
      // (30 km selected to allow misses of smaller size when rte_pos2 is at
      // surface level, but surface interference never OK if rte_pos above TOA)
      if( ppath_what_background(ppt) == 2  &&  ip == ppt.np-1  
                                                           &&  l12-lip > 30e3 )
        { 
          za_new       = rte_los[0] - 1; 
          za_upp_limit = rte_los[0]; 
        }

      // Ppath OK
      else
        {
          // Estimate ppath at the distance of l12, and calculate size
          // of "miss" (measured in diffference in geometric angles)
          Vector  los;
          Numeric dza;          
          if( atmosphere_dim <= 2 )
            { 
              // Convert pos and los for point ip to cartesian coordinates
              Numeric xip, zip;
              poslos2cart( xip, zip, dxip, dzip, ppt.r[ip], ppt.pos(ip,1), 
                                                              ppt.los(ip,0) );
              // Find where the extension from point ip crosses the l12 
              // sphere: point c
              Numeric latc;
              line_circle_intersect( xc, zc, xip, zip, dxip, dzip, x1, z1, l12);
              cart2pol( rc, latc, xc, zc, ppt.pos(ip,1), ppt.los(ip,0) );
              posc[1] = latc;
              posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                                    lat_grid, lon_grid, posc );
              rte_losGeometricFromRtePosToRtePos2( los, atmosphere_dim, 
                  lat_grid, lon_grid, refellipsoid, rte_pos, posc, verbosity );
              dza = los[0] - rte_los_geom[0];
            }
          else
            { 
              // Convert pos and los for point ip to cartesian coordinates
              Numeric xip, yip, zip;
              poslos2cart( xip, yip, zip, dxip, dyip, dzip, ppt.r[ip], 
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) );
              // Find where the extension from point ip crosses the l12 
              // sphere: point c
              Numeric latc, lonc;
              line_sphere_intersect( xc, yc, zc, xip, yip, zip, 
                                           dxip, dyip, dzip, x1, y1, z1, l12 );
              cart2sph( rc, latc, lonc, xc, yc, zc, ppt.pos(ip,1), 
                           ppt.pos(ip,2), ppt.los(ip,0), ppt.los(ip,1) );
              posc[1] = latc;
              posc[2] = lonc;
              posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                                    lat_grid, lon_grid, posc );
            }
          //
          rte_losGeometricFromRtePosToRtePos2( los, atmosphere_dim, 
                  lat_grid, lon_grid, refellipsoid, rte_pos, posc, verbosity );
          //
          dza = los[0] - rte_los_geom[0];

          // Update bookkeeping variables
          it++;
          t_za[it]  = rte_los[0];
          t_dza[it] = dza; 
          cout << "[\n";
          for( Index k=0; k<=it; k++ )
            {
              cout << fixed << setprecision(7) << t_za[k] << " " 
                   << t_dza[k] << endl;
            }
          cout << "];\n";
          //
          if( dza > 0  &&  rte_los[0] < za_upp_limit )
            { za_upp_limit = rte_los[0]; }
          else if( dza < 0  &&  rte_los[0] > za_low_limit )
            { za_low_limit = rte_los[0]; }
          cout << za_low_limit << " - " << za_upp_limit << endl;

          // Ready ?
          if( abs(dza) <= za_accuracy )
            { ready = true; }
          else if( za_upp_limit - za_low_limit <= za_accuracy/10 )
            { 
              failed = true; 
              cout << "Zenith angle search range closed !!!\n";
              break; 
            }
          // Catch non-convergence (just for extra safety, za-range should be
          // closed quicker than this)
          ntries += 1;
          if( ntries >= 100 )
            { 
              failed = true; 
              cout << "Too many iterations !!!\n";
              break; 
            }

          // Estimate new angle
          if( it < 1 )
            { za_new = rte_los[0] - dza; }
          else
            {
              // Estimate new angle by linear regression over some of the
              // last calculations
              const Index nfit = min( it+1, (Index)3 );
              const Index i0 = it - nfit + 1;
              Vector p;
              linreg( p, t_za[Range(i0,nfit)], t_dza[Range(i0,nfit)] );
              za_new = -p[0]/p[1];
            }
          //
          if( atmosphere_dim == 3 )
            { daa = los[1] - rte_los_geom[1]; }
        }

      // Update rte_losif not ready
      if( !ready && !failed )
        {
          // Update rte_los. Use bisection of za_new is basically
          // identical to old angle, or is outside lower or upper
          // limit. Otherwise use reult of linear reg.
          if( isinf(za_new)                                 ||  
              isnan(za_new)                                 ||  
              abs(za_new-rte_los[0]) < 0.99*za_accuracy ||
              za_new <= za_low_limit                        || 
              za_new >= za_upp_limit                         )
            { 
              rte_los[0] = (za_low_limit+za_upp_limit)/2;  
            }
          else
            { 
              rte_los[0] = za_new;
              if( atmosphere_dim == 3 )
                { rte_los[1] -= daa; }
            }
        }
      cout << "New rte_los: " << rte_los << endl;
    } // while
  //--------------------------------------------------------------------------

  if( failed )
    {
      throw runtime_error( "Some failure ...." );
    }


  // Create final ppath. 
  // If ground intersection: Set to length 1 and ground background,  
  // to flag non-OK path
  // Otherwise: Fill path and set background to transmitter

  if( ground )
    { 
      ppath_init_structure( ppath, atmosphere_dim, 1 ); 
      ppath_set_background( ppath, 2 );
    }

  else
    {
      // Distance between point ip of ppt and posc
      Numeric ll;
      if( atmosphere_dim <= 2 )
        { distance2D( ll, rc, posc[1], ppt.r[ip], ppt.pos(ip,1) ); }
      else
        { distance3D( ll, rc, posc[1], posc[2], ppt.r[ip], ppt.pos(ip,1), 
                                                           ppt.pos(ip,2) ); }

      // Last point of ppt closest to rte_pos2. No point to add, maybe
      // calculate start_lstep and start_los:
      if( ip == ppt.np-1 )
        { 
          ppath_init_structure( ppath, atmosphere_dim, ppt.np );
          ppath_copy( ppath, ppt, -1 );
          if( ppath_what_background( ppath ) == 1 )
            { 
              ppath.start_lstep= ll; 
              Numeric d1, d2=0, d3;
          if( atmosphere_dim <= 2 )
            { cart2poslos( d1, d3, ppath.start_los[0], xc, zc, dxip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.los(ip,0) ); }
          else
            { cart2poslos( d1, d2, d3, ppath.start_los[0], ppath.start_los[1],
                           xc, yc, zc, dxip, dyip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) ); }
            }
        }
      // rte_pos2 inside the atmosphere (posc entered as end point) 
      else
        {
          ppath_init_structure( ppath, atmosphere_dim, ip+2 );
          ppath_copy( ppath, ppt, ip+1 );
          //
          const Index i = ip+1;
          if( atmosphere_dim <= 2 )
            { cart2poslos( ppath.r[i], ppath.pos(i,1), ppath.los(i,0), xc, zc, 
                           dxip, dzip, ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.los(ip,0) ); }
          else
            { cart2poslos( ppath.r[i], ppath.pos(i,1), ppath.pos(i,2), 
                           ppath.los(i,0), ppath.los(i,1), xc, yc, zc, 
                           dxip, dyip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) ); }
          //
          ppath.pos(i,joker) = posc;
          ppath.lstep[i-1]   = ll;
          ppath.start_los    = ppath.los(i,joker);
          
          // n by linear interpolation
          assert( ll < ppt.lstep[i-1] );
          const Numeric w = ll/ppt.lstep[i-1];
          ppath.nreal[i]  = (1-w)*ppt.nreal[i-1]  + w*ppt.nreal[i];
          ppath.ngroup[i] = (1-w)*ppt.ngroup[i-1] + w*ppt.ngroup[i];

          // Grid positions
          GridPos gp_lat, gp_lon;
          rte_pos2gridpos( ppath.gp_p[i], gp_lat, gp_lon, 
                           atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                           ppath.pos(i,Range(0,atmosphere_dim)) );
          if( atmosphere_dim >= 2 )
            { 
              gridpos_copy( ppath.gp_lat[i], gp_lat );
              if( atmosphere_dim == 3 )
                { gridpos_copy( ppath.gp_lon[i], gp_lon ); }
            }
        }

      // Common stuff
      ppath_set_background( ppath, 9 );
      ppath.start_pos = rte_pos2;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ppathFromRtePos2V2Old(      
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
    const Tensor3&        edensity_field,
    const Vector&         f_grid,
    const Vector&         refellipsoid,
    const Matrix&         z_surface,
    const Vector&         rte_pos,
    const Vector&         rte_pos2,
    const Numeric&        ppath_lraytrace,
    const Verbosity&      verbosity )
{
  //--- Check input -----------------------------------------------------------
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_pos( atmosphere_dim, rte_pos2, true );
  //---------------------------------------------------------------------------

 
  // Radius of rte_pos and rte_pos2
  const Numeric r1 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                              lon_grid, rte_pos ) + rte_pos[0];
  const Numeric r2 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                            lon_grid, rte_pos2 ) + rte_pos2[0];

  // Geometric distance between rte_pos and rte_pos2, effective 2D-lat for
  // rte_pos and and Cartesian coordinates of rte_pos:
  Numeric l12, lat1=0, x1, y1=0, z1;
  if( atmosphere_dim <= 2 )
    { 
      if( atmosphere_dim == 2 )
        { lat1 = rte_pos[1]; }
      distance2D( l12, r1, lat1, r2, rte_pos2[1] ); 
      pol2cart( x1, z1, r1, lat1 );
    } 
  else
    { 
      distance3D( l12, r1, rte_pos[1],  rte_pos[2], 
                       r2, rte_pos2[1], rte_pos2[2] ); 
      sph2cart( x1, y1, z1, r1, rte_pos[1], rte_pos[2] );
    }

  // Geometric LOS from rte_pos to rte_pos2
  Vector rte_los_geom;
  rte_losGeometricFromRtePosToRtePos2( rte_los_geom, atmosphere_dim, lat_grid, 
                        lon_grid, refellipsoid, rte_pos, rte_pos2, verbosity );

  // Add handling of 2D with negative za !!!!!!

  // Define remaining variables used in the while-loop below
  //
  // Basic bookkeeping variables
  Numeric za_upp_limit = 180;
  Numeric za_low_limit = 0;
  //
  // First try rte_los is taken from the geometric path
  Vector rte_los_try = rte_los_geom;
  //
  // Various variables associated with the ppath, and the point of the path
  // closest to the transmitter
  Ppath   ppt;                // "Test ppath"
  Index   ip  = -999;         // Index of closest ppath point
  Numeric dxip, dyip=0, dzip; // Cartesian LOS of the closest ppath point
  //
  // Data for the intersection of the l12-sphere
  Vector  posc(max( Index(2),atmosphere_dim));
  Numeric rc, xc, yc=0, zc;   

  CREATE_OUT2;
  CREATE_OUT3;

  Vector t_za(99), t_dza(99);
  Index  it = -1;

  const Numeric za_accuracy = 1e-5;
  const Numeric nlfac = 2;

  // Keep trying until ready or ground intersetion determined
  //
  bool  ready  = false;
  bool  ground = false;
  bool  failed = false;
  Index ntries = 0;
  //
  while( !ready && !ground )
    {
      out2 << "    Trying rte_los: [" << rte_los_try << "]\n";

      // Path for present rte_los (no cloudbox!)
      ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
                  lon_grid, t_field, z_field, vmr_field, edensity_field,
                  f_grid, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
                  rte_pos, rte_los_try, ppath_lraytrace, 0, verbosity );

      // Find the point closest to rte_pos2, on the side towards rte_pos. 
      // We do this by looking at the distance to rte_pos, that should be 
      // as close to l12 as possible, but not exceed it.
      Numeric lip = 99e99;
      ip = ppt.np;
      //
      while( lip >= l12  &&  ip > 0 )
        {
          ip--;
          if( atmosphere_dim <= 2 )
            { distance2D( lip, r1, lat1, ppt.r[ip], ppt.pos(ip,1) ); }
          else
            { distance3D( lip, r1, rte_pos[1], rte_pos[2], 
                          ppt.r[ip], ppt.pos(ip,1), ppt.pos(ip,2) ); }
        }

      // Surface intersection:
      // Not OK if the ground position is too far from rte_pos2. 
      // (30 km selected to allow misses of smaller size when rte_pos2 is at
      // surface level, but surface interference never OK if rte_pos above TOA)
      if( ppath_what_background(ppt) == 2  &&  ip == ppt.np-1  
                                                           &&  l12-lip > 30e3 )
        { 
          if( rte_los_try[0] < za_upp_limit )
            { za_upp_limit = rte_los_try[0]; }

          // New za by bisection or a fixed step
          const Numeric ground_dza = 1;
          const Numeric za_unc     = rte_los_try[0] - za_low_limit;
          if( za_unc > 2*ground_dza )
            { rte_los_try[0] -= ground_dza; }
          else
            {
              if( za_unc <= za_accuracy/nlfac )
                { 
                  ground = true; 
                  cout << "Ground intersection confirmed !!!\n";
                }
              else
                { rte_los_try[0] -= za_unc/2.0; }
            }
        }

      // Ppath OK
      else
        {
          // Estimate ppath at the distance of l12, and calculate size
          // of "miss" (measured in diffference in geometric angles)
          Vector  los;
          Numeric dza;          
          if( atmosphere_dim <= 2 )
            { 
              // Convert pos and los for point ip to cartesian coordinates
              Numeric xip, zip;
              poslos2cart( xip, zip, dxip, dzip, ppt.r[ip], ppt.pos(ip,1), 
                                                              ppt.los(ip,0) );
              // Find where the extension from point ip crosses the l12 
              // sphere: point c
              Numeric latc;
              line_circle_intersect( xc, zc, xip, zip, dxip, dzip, x1, z1, l12);
              cart2pol( rc, latc, xc, zc, ppt.pos(ip,1), ppt.los(ip,0) );
              posc[1] = latc;
              posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                                    lat_grid, lon_grid, posc );
              rte_losGeometricFromRtePosToRtePos2( los, atmosphere_dim, 
                  lat_grid, lon_grid, refellipsoid, rte_pos, posc, verbosity );
              dza = los[0] - rte_los_geom[0];
            }
          else
            { 
              // Convert pos and los for point ip to cartesian coordinates
              Numeric xip, yip, zip;
              poslos2cart( xip, yip, zip, dxip, dyip, dzip, ppt.r[ip], 
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) );
              // Find where the extension from point ip crosses the l12 
              // sphere: point c
              Numeric latc, lonc;
              line_sphere_intersect( xc, yc, zc, xip, yip, zip, 
                                           dxip, dyip, dzip, x1, y1, z1, l12 );
              cart2sph( rc, latc, lonc, xc, yc, zc, ppt.pos(ip,1), 
                           ppt.pos(ip,2), ppt.los(ip,0), ppt.los(ip,1) );
              posc[1] = latc;
              posc[2] = lonc;
              posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                                    lat_grid, lon_grid, posc );
            }
          //
          rte_losGeometricFromRtePosToRtePos2( los, atmosphere_dim, 
                  lat_grid, lon_grid, refellipsoid, rte_pos, posc, verbosity );


          // Ready or update rte_los_try
          //
          dza = los[0] - rte_los_geom[0];
          //
          it++;
          t_za[it]  = rte_los_try[0];
          t_dza[it] = dza; 
          cout << "[\n";
          for( Index k=0; k<=it; k++ )
            {
              cout << fixed << setprecision(7) << t_za[k] << " " 
                   << t_dza[k] << endl;
            }
          cout << "];\n";
          //
          if( abs(dza) <= za_accuracy )
            { ready = true; }
          else
            {
              if( dza > 0  &&  rte_los_try[0] < za_upp_limit )
                { za_upp_limit = rte_los_try[0]; }
              else if( dza < 0  &&  rte_los_try[0] > za_low_limit )
                { za_low_limit = rte_los_try[0]; }

              cout << za_low_limit << " - " << za_upp_limit << endl;

              if( za_upp_limit - za_low_limit <= za_accuracy/nlfac )
                { 
                  failed = true; 
                  cout << "Zenith angle search range closed !!!\n";
                  break; 
                }


              //it++;
              //t_za[it]  = rte_los_try[0];
              //t_dza[it] = dza; 

              if( it < 1 )
                {
                  if( dza > 0 )
                    {
                      if( rte_los_try[0]-dza > za_low_limit )
                        { rte_los_try[0] -= dza; }
                      else
                        { rte_los_try[0] = (rte_los_try[0]+za_low_limit)/2; }
                    }
                  else
                    {
                      if( rte_los_try[0]-dza < za_upp_limit )
                        { rte_los_try[0] -= dza; }
                      else
                        { rte_los_try[0] = (rte_los_try[0]+za_upp_limit)/2; }
                    }
                  if( atmosphere_dim == 3 )
                    { rte_los_try[1] -= los[1] - rte_los_geom[1]; }
                }
              else
                {
                  // Estimate new angle by linear regression over some of the
                  // last calculations
                  const Index nfit = min( it+1, (Index)3 );
                  const Index i0 = it - nfit + 1;
                  Vector p;
                  linreg( p, t_za[Range(i0,nfit)], t_dza[Range(i0,nfit)] );
                  const Numeric za_new = -p[0]/p[1];

                  // Update rte_los_try. Use bisection of za_new is basically
                  // identical to old angle, or is outside lower or upper
                  // limit. Otherwise use reult of linear reg.
                  if( isinf(za_new) ||  isnan(za_new)  ||   
                      abs(za_new-rte_los_try[0]) < 0.99*za_accuracy )
                    { rte_los_try[0] = (za_low_limit+za_upp_limit)/2; }
                  else if( za_new <= za_low_limit )
                    { rte_los_try[0] = (rte_los_try[0]+za_low_limit)/2; }
                  else if( za_new >= za_upp_limit )
                    { rte_los_try[0] = (rte_los_try[0]+za_upp_limit)/2; }
                  else
                    { 
                      rte_los_try[0] = za_new;
                      if( atmosphere_dim == 3 )
                        { rte_los_try[1] -= los[1] - rte_los_geom[1]; }
                    }
                }
              cout << "New rte_los:\n" << rte_los_try << endl;
            }
        }  // else

      // Catch non-convergence (just for extra safety, za-range should be
      // closed quicker than this)
      ntries += 1;
      if( ntries >= 100 )
        { 
          failed = true; 
          cout << "Too many iterations !!!\n";
          break; 
        }
    } // while
  //--------------------------------------------------------------------------

  if( failed )
    {
      throw runtime_error( "Some failure ...." );
    }


  // Create final ppath. 
  // If ground intersection: Set to length 1 and ground background,  
  // to flag non-OK path
  // Otherwise: Fill path and set background to transmitter

  if( ground )
    { 
      ppath_init_structure( ppath, atmosphere_dim, 1 ); 
      ppath_set_background( ppath, 2 );
    }

  else
    {
      // Distance between point ip of ppt and posc
      Numeric ll;
      if( atmosphere_dim <= 2 )
        { distance2D( ll, rc, posc[1], ppt.r[ip], ppt.pos(ip,1) ); }
      else
        { distance3D( ll, rc, posc[1], posc[2], ppt.r[ip], ppt.pos(ip,1), 
                                                           ppt.pos(ip,2) ); }

      // Last point of ppt closest to rte_pos2. No point to add, maybe
      // calculate start_lstep and start_los:
      if( ip == ppt.np-1 )
        { 
          ppath_init_structure( ppath, atmosphere_dim, ppt.np );
          ppath_copy( ppath, ppt, -1 );
          if( ppath_what_background( ppath ) == 1 )
            { 
              ppath.start_lstep= ll; 
              Numeric d1, d2=0, d3;
          if( atmosphere_dim <= 2 )
            { cart2poslos( d1, d3, ppath.start_los[0], xc, zc, dxip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.los(ip,0) ); }
          else
            { cart2poslos( d1, d2, d3, ppath.start_los[0], ppath.start_los[1],
                           xc, yc, zc, dxip, dyip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) ); }
            }
        }
      // rte_pos2 inside the atmosphere (posc entered as end point) 
      else
        {
          ppath_init_structure( ppath, atmosphere_dim, ip+2 );
          ppath_copy( ppath, ppt, ip+1 );
          //
          const Index i = ip+1;
          if( atmosphere_dim <= 2 )
            { cart2poslos( ppath.r[i], ppath.pos(i,1), ppath.los(i,0), xc, zc, 
                           dxip, dzip, ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.los(ip,0) ); }
          else
            { cart2poslos( ppath.r[i], ppath.pos(i,1), ppath.pos(i,2), 
                           ppath.los(i,0), ppath.los(i,1), xc, yc, zc, 
                           dxip, dyip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) ); }
          //
          ppath.pos(i,joker) = posc;
          ppath.lstep[i-1]   = ll;
          ppath.start_los    = ppath.los(i,joker);
          
          // n by linear interpolation
          assert( ll < ppt.lstep[i-1] );
          const Numeric w = ll/ppt.lstep[i-1];
          ppath.nreal[i]  = (1-w)*ppt.nreal[i-1]  + w*ppt.nreal[i];
          ppath.ngroup[i] = (1-w)*ppt.ngroup[i-1] + w*ppt.ngroup[i];

          // Grid positions
          GridPos gp_lat, gp_lon;
          rte_pos2gridpos( ppath.gp_p[i], gp_lat, gp_lon, 
                           atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                           ppath.pos(i,Range(0,atmosphere_dim)) );
          if( atmosphere_dim >= 2 )
            { 
              gridpos_copy( ppath.gp_lat[i], gp_lat );
              if( atmosphere_dim == 3 )
                { gridpos_copy( ppath.gp_lon[i], gp_lon ); }
            }
        }

      // Common stuff
      ppath_set_background( ppath, 9 );
      ppath.start_pos = rte_pos2;
    }

  rte_los = rte_los_try;
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
    const Tensor3&        edensity_field,
    const Vector&         f_grid,
    const Vector&         refellipsoid,
    const Matrix&         z_surface,
    const Vector&         rte_pos,
    const Vector&         rte_pos2,
    const Numeric&        ppath_lraytrace,
    const Numeric&        dl_target,
    const Numeric&        dl_ok,
    const Verbosity&      verbosity )
{
  //--- Check input -----------------------------------------------------------
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_pos( atmosphere_dim, rte_pos2, true );
  //---------------------------------------------------------------------------

  CREATE_OUT2;
  CREATE_OUT3;

  // Radius of rte_pos and rte_pos2
  const Numeric r1 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                              lon_grid, rte_pos ) + rte_pos[0];
  const Numeric r2 = pos2refell_r( atmosphere_dim, refellipsoid, lat_grid, 
                                            lon_grid, rte_pos2 ) + rte_pos2[0];
  
  // Geometric LOS from rte_pos to rte_pos2
  Vector rte_los_geom;
  rte_losGeometricFromRtePosToRtePos2( rte_los_geom, atmosphere_dim, lat_grid, 
                        lon_grid, refellipsoid, rte_pos, rte_pos2, verbosity );

  // Geometric distance between rte_pos and rte_pos2, effective 2D-lat for
  // rte_pos and and Cartesian coordinates of rte_pos:
  Numeric l12, lat1=0, x1, y1=0, z1;
  if( atmosphere_dim <= 2 )
    { 
      if( atmosphere_dim == 2 )
        { lat1 = rte_pos[1]; }
      distance2D( l12, r1, lat1, r2, rte_pos2[1] ); 
      pol2cart( x1, z1, r1, lat1 );
    } 
  else
    { 
      distance3D( l12, r1, rte_pos[1],  rte_pos[2], 
                       r2, rte_pos2[1], rte_pos2[2] ); 
      sph2cart( x1, y1, z1, r1, rte_pos[1], rte_pos[2] );
    }

  // A matrix with "hit data", to keep track of most close tries so far.
  // Row 0 holds data for the closest calculation where end point is at a too
  // high zenith angle. Row 1, the same for too low zenith angle. The columns,
  // 0: Valid data flag. Below 0, not data yet. Above 0, valid data.
  // 1: Size of zenith angle deviation for end point
  // 2: Zenith angle of this calculation
  // 3: Estimate of path length error
  // 4. Spatial error.
  // 5. Azimuth angle.
  Matrix  H(2,6);
  H(joker,0) = -1; 
  H(0,1) = 99e99;    // A value to indicate infinite error
  H(1,1) = -99e99;   // A value to indicate infinite error
  H(0,3) = 99e99;    // A value to indicate infinite error
  H(1,3) = 99e99;    // A value to indicate infinite error

  // Iterate until convergence for rte_los
  rte_los = rte_los_geom;  // We start with geomtrical angles
  Index   ready = false;   // True when convergence criteria fulfilled
  Index   counter=0;       // Number of tries
  Numeric dl=1e9;          // Path length error
  Numeric ds;              // Spatial error
  Index   surface_hits=0;  // Number of tries ending up in the surface
  Ppath   ppt, ppt_old;    // Present and old test ppath
  // The point of ppath closest to rte_pos2 (on "l12 sphere"), and 
  // its index, radius and cartesian coordinates
  Vector  posc( max( Index(2), atmosphere_dim ) );
  Index   ip=-999;
  Numeric rc, xc, yc=0, zc, dxip, dyip=0, dzip;
  // The last variables must be remembered for old test path
  Index   ip_old = -1;
  Numeric rc_old = -999, xc_old   = -999, yc_old   = -999, zc_old   = -999;
  Numeric dl_old = -999, dxip_old = -999, dyip_old = -999, dzip_old = -999;
  //
  while( !ready )
    {
      out3 << "    Trying rte_los: [" << rte_los << "]\n";

      // Save old ppt stuff
      if( counter )
        {
          ppath_init_structure( ppt_old, atmosphere_dim, ppt.np );
          ppath_copy( ppt_old, ppt, ppt.np );
          ip_old   = ip; 
          dl_old   = dl;
          rc_old   = rc; 
          xc_old   = xc; 
          yc_old   = yc; 
          zc_old   = zc;
          dxip_old = dxip;   
          dyip_old = dyip;   
          dzip_old = dzip;
        }

      // Path for present rte_los (no cloudbox!)
      ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
                  lon_grid, t_field, z_field, vmr_field, edensity_field,
                  f_grid, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
                  rte_pos, rte_los, ppath_lraytrace, 0, verbosity );

      // Find the point closest to rte_pos2, on the side towards rte_pos. 
      // We do this by looking at the distance to rte_pos, that should be 
      // as close to l12 as possible, but not exceed it.
      Numeric lip = 99e99;
      ip = ppt.np;
      //
      while( lip >= l12  &&  ip > 0 )
        {
          ip--;
          
          // Distance between rte_pos and ppath point
          if( atmosphere_dim <= 2 )
            { distance2D( lip, r1, lat1, ppt.r[ip], ppt.pos(ip,1) ); }
          else
            { distance3D( lip, r1, rte_pos[1], rte_pos[2], 
                          ppt.r[ip], ppt.pos(ip,1), ppt.pos(ip,2) ); }
        }

      // Surface intersections too far from rte_pos2 not OK 
      // (50 km selected to allow misses of smaller size when rte_pos2 is at
      // surface level, but surface interference never OK if rte_pos above TOA)
      if( ppath_what_background(ppt) == 2  &&  ip == ppt.np-1  
                                                           &&  l12-lip > 50e3 )
        { 
          //cout << "Surface hit\n";
          // Stop if "too many" surface hits
          surface_hits++;
          if( surface_hits == 20 )
            { 
              ostringstream os;
              os << "Repeated intersections with the surface caused failure\n"
                 << "in the attempt to determine propagation path.";
              throw runtime_error( os.str() );
            }

          // Set H, with a relative large error
          H(0,0)=1;  H(0,1)= 1/sqrt(surface_hits) ; 
          H(0,2)=rte_los[0]; H(0,3)=99e3; H(0,4)=99e3; 
          if( atmosphere_dim == 3 ) { H(0,5)=rte_los[1]; }

          // If an OK lower angle exist, take middle point.
          // Otherwise, make a substantial jump upwards
          if( H(1,0) > 0 )
            { rte_los[0] =  0.5*rte_los[0] + 0.5*H(1,2); }
          else
            { rte_los[0] -= 0.5; }
          //cout << "-----\n" << H << "\n-----\n" << endl;
        }

      // No surface intersection:
      else
        {
          // Estimate ppath at the distance of l12, and the distance from
          // that point to rte_pos2 (ds)
          if( atmosphere_dim <= 2 )
            { 
              // Convert pos and los for point ip to cartesian coordinates
              Numeric xip, zip;
              poslos2cart( xip, zip, dxip, dzip, ppt.r[ip], ppt.pos(ip,1), 
                                                              ppt.los(ip,0) );
              // Find where the extension from point ip crosses the l12 
              // sphere: point c
              Numeric latc;
              line_circle_intersect( xc, zc, xip, zip, dxip, dzip, x1, z1, l12);
              cart2pol( rc, latc, xc, zc, ppt.pos(ip,1), ppt.los(ip,0) );
              posc[1] = latc;
              posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                                    lat_grid, lon_grid, posc );
              distance2D( ds, rc, latc, r2, rte_pos2[1] );
              dl = sqrt( l12*l12 + ds*ds ) - l12;
            }
          else
            { 
              // Convert pos and los for point ip to cartesian coordinates
              Numeric xip, yip, zip;
              poslos2cart( xip, yip, zip, dxip, dyip, dzip, ppt.r[ip], 
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) );
              // Find where the extension from point ip crosses the l12 
              // sphere: point c
              Numeric latc, lonc;
              line_sphere_intersect( xc, yc, zc, xip, yip, zip, 
                                           dxip, dyip, dzip, x1, y1, z1, l12 );
              cart2sph( rc, latc, lonc, xc, yc, zc, ppt.pos(ip,1), 
                           ppt.pos(ip,2), ppt.los(ip,0), ppt.los(ip,1) );
              posc[1] = latc;
              posc[2] = lonc;
              posc[0] = rc - pos2refell_r( atmosphere_dim, refellipsoid, 
                                                    lat_grid, lon_grid, posc );
              distance3D( ds, rc, latc, lonc, r2, rte_pos2[1], rte_pos2[2] );
              dl = sqrt( l12*l12 + ds*ds ) - l12;
            }
          //cout << "-----\n" << H << "\n-----\n" << endl;
          

          // Converged?
          if( dl <= dl_target )
            { ready = true;  } 
          
          // Otherwise, calculate new rte_los
          else
            {
              Vector los;
              rte_losGeometricFromRtePosToRtePos2( los, atmosphere_dim, 
                  lat_grid, lon_grid, refellipsoid, rte_pos, posc, verbosity );
              Numeric dza = los[0] - rte_los_geom[0];
              //cout << rte_los[0] - rte_los_geom[0] << " " << dza << endl;

              // Any improvement compared to existing results in H
              //
              bool improving = true;
              //
              if( dza > 0   &&  dza < H(0,1) )
                { 
                  H(0,0)=1;  H(0,1)=dza; H(0,2)=rte_los[0]; 
                  H(0,3)=dl; H(0,4)=ds; 
                  if( atmosphere_dim == 3 ) { H(0,5) = rte_los[1]; }
                }
              else if( dza < 0   &&  dza > H(1,1) )
                { 
                  H(1,0)=1;  H(1,1)=dza; H(1,2)=rte_los[0]; 
                  H(1,3)=dl; H(1,4)=ds; 
                  if( atmosphere_dim == 3 ) { H(1,5) = rte_los[1]; }
                } 
              else // Ending up here means "jitter" in ppath calculation
                { improving = false; }

              if( improving )
                {
                  // Select new rte_los. If both rows of H filled, use
                  // bisection. Otherwise use dza.
                  if( H(0,0) > 0   &&   H(1,0) > 0 )
                    { rte_los[0] = H(1,2) - (H(0,2)-H(1,2)) * H(1,1) / 
                                                             (H(0,1)-H(1,1)); }
                  else
                    { rte_los[0] -= dza; }

                  // For azimuth angle, we just apply the "miss angle"
                  if( atmosphere_dim == 3 )
                    { rte_los[1] -= los[1] - rte_los_geom[1]; }
                }
              else
                {
                  // OK or fail?
                  if( dl_old <= dl_ok )
                    { 
                      // Switch back to old ppt and rc, and we are ready
                      ppath_init_structure( ppt, atmosphere_dim, ppt_old.np );
                      ppath_copy( ppt, ppt_old, ppt_old.np );
                      ip    = ip_old; 
                      dl    = dl_old;
                      rc    = rc_old; 
                      xc    = xc_old; 
                      yc    = yc_old; 
                      zc    = zc_old;
                      dxip  = dxip_old;   
                      dyip  = dyip_old;   
                      dzip  = dzip_old;
                      ready = true;  
                    } 
                  else
                    {
                      if( H(0,3) < H(1,3) )  { dl = H(0,3);   ds = H(0,4); }
                      else                   { dl = H(1,3);   ds = H(1,4); }
                      ostringstream os;
                      os << "The smallest path length error achieved is: " 
                         << dl_old
                         << " m\nwhich is above the set accuracy limit of  : " 
                         << dl_ok << " m.\nGeographically, the closest "
                         << "match with *rte_pos2* is " << ds << " m.\n"
                         << "This can potentially be improved by decreasing the"
                         << " value of\n*ppath_lraytrace*, and maybe also the "
                         << "one of *ppath_lmax*.\n";
                      throw runtime_error( os.str() );
                    }
                }
            }
        }

      // Brake if not finding a solution
      counter++;
      if( counter >= 50 )
        throw runtime_error( "Sorry, but the algorithm is not converging. "
                             "Don't know what to do!" );
    }

  out2 << "  Length error: " << dl << " m.\n";
  out2 << "  Spatial miss: " << ds << " m.\n";

  // Create final ppath. BAckground set to "transmitter".

  // Distance between point ip of ppt and posc
  Numeric ll;
  if( atmosphere_dim <= 2 )
    { distance2D( ll, rc, posc[1], ppt.r[ip], ppt.pos(ip,1) ); }
  else
    { distance3D( ll, rc, posc[1], posc[2], ppt.r[ip], ppt.pos(ip,1), 
                                                       ppt.pos(ip,2) ); }

  // Last point of ppt closest to rte_pos2. No point to add, maybe calculate
  // start_lstep and start_los: 
  if( ip == ppt.np-1 )
    { 
      ppath_init_structure( ppath, atmosphere_dim, ppt.np );
      ppath_copy( ppath, ppt, -1 );
      if( ppath_what_background( ppath ) == 1 )
        { 
          ppath.start_lstep= ll; 
          Numeric d1, d2=0, d3;
          if( atmosphere_dim <= 2 )
            { cart2poslos( d1, d3, ppath.start_los[0], xc, zc, dxip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.los(ip,0) ); }
          else
            { cart2poslos( d1, d2, d3, ppath.start_los[0], ppath.start_los[1],
                           xc, yc, zc, dxip, dyip, dzip, 
                           ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                           ppt.pos(ip,1), ppt.pos(ip,2), 
                           ppt.los(ip,0), ppt.los(ip,1) ); }
        }
    }
  // rte_pos2 inside the atmosphere (posc entered as end point) 
  else
    {
      ppath_init_structure( ppath, atmosphere_dim, ip+2 );
      ppath_copy( ppath, ppt, ip+1 );
      //
      const Index i = ip+1;
      if( atmosphere_dim <= 2 )
        { cart2poslos( ppath.r[i], ppath.pos(i,1), ppath.los(i,0), xc, zc, 
                       dxip, dzip, ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                       ppt.pos(ip,1), ppt.los(ip,0) ); }
      else
        { cart2poslos( ppath.r[i], ppath.pos(i,1), ppath.pos(i,2), 
                       ppath.los(i,0), ppath.los(i,1), xc, yc, zc, 
                       dxip, dyip, dzip, ppt.r[ip]*sin(DEG2RAD*ppt.los(ip,0)),
                       ppt.pos(ip,1), ppt.pos(ip,2), 
                       ppt.los(ip,0), ppt.los(ip,1) ); }
      //
      ppath.pos(i,joker) = posc;
      ppath.lstep[i-1]   = ll;
      ppath.start_los    = ppath.los(i,joker);

      // n by linear interpolation
      assert( ll < ppt.lstep[i-1] );
      const Numeric w = ll/ppt.lstep[i-1];
      ppath.nreal[i]  = (1-w)*ppt.nreal[i-1]  + w*ppt.nreal[i];
      ppath.ngroup[i] = (1-w)*ppt.ngroup[i-1] + w*ppt.ngroup[i];

      // Grid positions
      GridPos gp_lat, gp_lon;
      rte_pos2gridpos( ppath.gp_p[i], gp_lat, gp_lon, 
                       atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                       ppath.pos(i,Range(0,atmosphere_dim)) );
      if( atmosphere_dim >= 2 )
        { 
          gridpos_copy( ppath.gp_lat[i], gp_lat );
          if( atmosphere_dim == 3 )
            { gridpos_copy( ppath.gp_lon[i], gp_lon ); }
        }
    }

  // Common stuff
  ppath_set_background( ppath, 9 );
  ppath.start_pos = rte_pos2;
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
    const Tensor3&        edensity_field,
    const Vector&         f_grid,
    const Vector&         refellipsoid,
    const Matrix&         z_surface,
    const Index&          cloudbox_on, 
    const ArrayOfIndex&   cloudbox_limits,
    const Vector&         rte_pos,
    const Vector&         rte_los,
    const Numeric&        ppath_lraytrace,
    const Verbosity&      verbosity)
{
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid, 
              lon_grid, t_field, z_field, vmr_field, edensity_field, f_grid, 
              refellipsoid, z_surface, cloudbox_on, cloudbox_limits, rte_pos, 
              rte_los, ppath_lraytrace, ppath_inside_cloudbox_do, verbosity );
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
      ppath_step.nreal[0]  = 1;
      ppath_step.ngroup[0] = 1;
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepRefractionBasic(
          Workspace&  ws,
          Ppath&      ppath_step,
    const Agenda&     refr_index_agenda,
    const Index&      atmosphere_dim,
    const Vector&     p_grid,
    const Vector&     lat_grid,
    const Vector&     lon_grid,
    const Tensor3&    z_field,
    const Tensor3&    t_field,
    const Tensor4&    vmr_field,
    const Tensor3&    edensity_field,
    const Vector&     refellipsoid,
    const Matrix&     z_surface,
    const Vector&     f_grid,
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
                              t_field, vmr_field, edensity_field, f_grid, 
                              refellipsoid, z_surface(0,0), ppath_lmax, 
                              refr_index_agenda, "linear_basic", 
                              ppath_lraytrace );
        }
      else if( atmosphere_dim == 2 )
        { 
          ppath_step_refr_2d( ws, ppath_step, p_grid, lat_grid, 
                              z_field(joker,joker,0), t_field, vmr_field, 
                              edensity_field, f_grid, refellipsoid, 
                              z_surface(joker,0), ppath_lmax, 
                              refr_index_agenda, "linear_basic", 
                              ppath_lraytrace ); 
        }
      else if( atmosphere_dim == 3 )
        { 
          ppath_step_refr_3d( ws, ppath_step, p_grid, lat_grid, lon_grid, 
                              z_field, t_field, vmr_field, edensity_field, 
                              f_grid, refellipsoid, z_surface, ppath_lmax, 
                              refr_index_agenda, "linear_basic", 
                              ppath_lraytrace ); 
        }
      else
        { throw runtime_error( "The atmospheric dimensionality must be 1-3." );}
    }

  else
    { 
      assert( ppath_step.np == 1 );
      if( atmosphere_dim == 1 )
        { get_refr_index_1d( ws, ppath_step.nreal[0], ppath_step.ngroup[0], 
                             refr_index_agenda, p_grid, refellipsoid, 
                             z_field(joker,0,0), t_field, vmr_field, 
                             edensity_field, f_grid, ppath_step.r[0] ); 
        }
      else if( atmosphere_dim == 2 )
        { get_refr_index_2d( ws, ppath_step.nreal[0], ppath_step.ngroup[0], 
                             refr_index_agenda, p_grid, lat_grid, refellipsoid, 
                             z_field(joker,joker,0), t_field, vmr_field, 
                             edensity_field, f_grid, ppath_step.r[0], 
                             ppath_step.pos(0,1) ); 
        }
      else
        { get_refr_index_3d( ws, ppath_step.nreal[0], ppath_step.ngroup[0], 
                             refr_index_agenda, p_grid, lat_grid, lon_grid, 
                             refellipsoid, z_field, t_field, vmr_field, 
                             edensity_field, f_grid, ppath_step.r[0], 
                             ppath_step.pos(0,1), ppath_step.pos(0,2) ); 
        }
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losSet(
          Vector&    rte_los,
    const Index&     atmosphere_dim,
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
  chk_rte_pos( atmosphere_dim, rte_pos );
  chk_rte_pos( atmosphere_dim, rte_pos2, true );

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
void rte_posSet(
          Vector&    rte_pos,
    const Index&     atmosphere_dim,
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
void rte_pos_losMoveToStartOfPpath(
          Vector&    rte_pos,
          Vector&    rte_los,
    const Index&     atmosphere_dim,
    const Ppath&     ppath,
    const Verbosity&)
{
  const Index np = ppath.np;

  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  if( np == 0 )
    throw runtime_error( "The input *ppath* is empty." );
  if( ppath.pos.nrows() != np )
    throw runtime_error( "Internal inconsistency in *ppath* (size of data " 
                         "does not match np)." );
  
  rte_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  if( atmosphere_dim < 3 )
    { rte_los = ppath.los(np-1,Range(0,1)); }
  else
    { rte_los = ppath.los(np-1,Range(0,2)); }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void TangentPointExtract(
          Vector&    tan_pos,
    const Ppath&     ppath,
    const Verbosity& )
{
  Index it;
  find_tanpoint( it, ppath );

  tan_pos.resize( ppath.pos.ncols() );

  if( it < 0 )
    {
      tan_pos = sqrt( -1 );  // = NaN
    }
  else
    {
      tan_pos[0] = ppath.pos(it,0);
      tan_pos[1] = ppath.pos(it,1);
      if( ppath.pos.ncols() == 3 )
        { tan_pos[2] = ppath.pos(it,2); }
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void TangentPointPrint(
    const Ppath&     ppath,
    const Index&     level,
    const Verbosity& verbosity)
{
  Index it;
  find_tanpoint( it, ppath );

  ostringstream os;

  if( it < 0 )
    {
      os << "Lowest altitude found at the end of the propagation path.\n"
         << "This indicates that the tangent point is either above the\n"
         << "top-of-the-atmosphere or below the planet's surface.";
    }
  else
    {
      os << "Tangent point position:\n-----------------------\n"
         << "     z = " << ppath.pos(it,0)/1e3 << " km\n" 
         << "   lat = " << ppath.pos(it,1) << " deg";
        if( ppath.pos.ncols() == 3 )
          os << "\n   lon: " << ppath.pos(it,2) << " deg";
    }

  CREATE_OUTS;
  SWITCH_OUTPUT (level, os.str ());  
}




/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZaRefr1D(
           Workspace&    ws,
           Vector&       za_vector,
     const Agenda&       refr_index_agenda,
     const Matrix&       sensor_pos,
     const Vector&       p_grid,
     const Tensor3&      t_field,
     const Tensor3&      z_field,
     const Tensor4&      vmr_field,
     const Tensor3&      edensity_field,
     const Vector&       refellipsoid,
     const Index&        atmosphere_dim,
     const Vector&       f_grid,
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
  Numeric   refr_index, refr_index_group;

  // Calculate refractive index for the tangential altitudes
  for( Index i=0; i<ztan_vector.nelem(); i++ ) 
    {
      get_refr_index_1d( ws, refr_index, refr_index_group, refr_index_agenda, 
                         p_grid, refellipsoid[0], z_field(joker,0,0), t_field, 
                         vmr_field, edensity_field, f_grid,
                         ztan_vector[i] + refellipsoid[0] );

    // Calculate zenith angle
    za_vector[i] = 180 - RAD2DEG* asin( refr_index * 
                                        (refellipsoid[0] + ztan_vector[i]) / 
                                        (refellipsoid[0] + sensor_pos(i,0)) );
  }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZa1D(
          Vector&       za_vector,
    const Matrix&       sensor_pos,
    const Vector&       refellipsoid,
    const Index&        atmosphere_dim,
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


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathWriteXMLPartial (//WS Input:
                           const String& file_format,
                           const Ppath&  ppath,
                           // WS Generic Input:
                           const String& f,
                           const Index&  file_index,
                           const Verbosity& verbosity)
{
    String filename = f;
    Ppath ppath_partial = ppath;
    ArrayOfGridPos empty_gp;
    //Vector empty_v;

    ppath_partial.gp_p = empty_gp;
    ppath_partial.gp_lat = empty_gp;
    ppath_partial.gp_lon = empty_gp;
    //ppath_partial.nreal = empty_v;
    //ppath_partial.ngroup = empty_v;

    if (file_index >= 0)
    {
        // Create default filename if empty
        filename_xml_with_index( filename, file_index, "ppath" );
    }

    WriteXML( file_format, ppath_partial, filename, "ppath", "", verbosity );
}



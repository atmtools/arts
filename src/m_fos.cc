/* Copyright (C) 2010
   Patrick Eriksson <patrick.eriksson@chalmers.se>
                            
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
  \file   m_fos.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2010-04-06 

  \brief  Workspace functions associated with the FOS scattering scheme.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "montecarlo.h"
#include "rte.h"
#include "special_interp.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;





/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void fos_yStandard(
        Workspace&                     ws,
        Tensor3&                       fos_y,
        Matrix&                        iy_error,
        Index&                         iy_error_type,
        Tensor3&                       iy_aux,
        ArrayOfTensor3&                diy_dx,
  const Vector&                        rte_pos,
  const Index&                         iy_aux_do,
  const Index&                         atmosphere_dim,
  const Vector&                        p_grid,
  const Vector&                        lat_grid,
  const Vector&                        lon_grid,
  const Tensor3&                       z_field,
  const Tensor3&                       t_field,
  const Tensor4&                       vmr_field,
  const Matrix&                        r_geoid,
  const Matrix&                        z_surface,
  const Index&                         cloudbox_on,
  const ArrayOfIndex&                  cloudbox_limits,
  const Index&                         stokes_dim,
  const Vector&                        f_grid,
  const Agenda&                        ppath_step_agenda,
  const Agenda&                        emission_agenda,
  const Agenda&                        abs_scalar_gas_agenda,
  const Agenda&                        iy_clearsky_agenda,
  const Tensor3&                       iy_transmission,
  const Tensor4&                       pnd_field,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                        opt_prop_gas_agenda,
  const Agenda&                        fos_y_agenda,
  const Matrix&                        fos_angles,
  const Index&                         fos_use_mean_scat_data,
  const Index&                         fos_n,
  const Index&                         fos_i )
{
  // Angles inside these ranges are considered to be equal for 1D and 2D
  const Numeric dza = 0.01;
  const Numeric daa = 1;

  const Index jacobian_do = 0;

  const Index   nfosa = fos_angles.nrows();
        Matrix  tmp;

  fos_y.resize(nfosa,f_grid.nelem(),stokes_dim);

  if( fos_i == fos_n-1 )
    {
      if( atmosphere_dim == 1 )
        { 
          for( Index ia=0; ia<nfosa; ia++ )
            { 
              // To use already calculated data we demand that difference
              // in za is < dza
              Index ihit = -1;
              //
              for( Index it=ia-1; it>=0 && ihit<0; it-- )
                {
                   if( fabs( fos_angles(ia,0) - fos_angles(it,0) ) < dza )
                     { ihit = it; }
                } 
              
              if( ihit >= 0 )
                { 
                  fos_y(ia,joker,joker) = fos_y(ihit,joker,joker); 
                }
              else
                {
                  iy_clearsky_agendaExecute( ws, tmp, iy_error, iy_error_type,
                                             iy_aux, diy_dx, 0, 
                                             rte_pos, fos_angles(ia,Range(0,1)),
                                             iy_transmission, 0, jacobian_do, 
                                             iy_aux_do, f_grid, p_grid, 
                                             lat_grid, lon_grid, t_field, 
                                             vmr_field, iy_clearsky_agenda );
                  fos_y(ia,joker,joker) = tmp;
                }
            }
          
        }
      else if( atmosphere_dim == 2 )
        { 
          Vector rte_los(1);

          for( Index ia=0; ia<nfosa; ia++ )
            { 
              // To use already calculated data we demand that difference
              // in za is < dza, and aa is mirrored with respect to the
              // orbit plane (aa(it) = -aa(ia)).
              Index ihit = -1;
              //
              for( Index it=ia-1; it>=0 && ihit<0; it-- )
                {
                   if( fabs( fos_angles(ia,0) - fos_angles(it,0) ) < dza )
                     { 
                       if( fabs( fos_angles(ia,1) + fos_angles(it,1) ) < daa )
                         { ihit = it; }
                     }
                } 
              
              if( ihit >= 0 )
                { 
                  fos_y(ia,joker,joker) = fos_y(ihit,joker,joker); 
                }
              else
                {
                  // LOS
                  if( fabs(fos_angles(ia,1)) <= 90 )
                    { rte_los[0] = fos_angles(ia,0); }
                  else
                    { rte_los[0] = -fos_angles(ia,0); }

                  // Create stretched latitude grid
                  //
                  Vector lat_stretched( lat_grid.nelem() );
                  //
                  // No strect needed for zenith, nadir and aa= 0 or +-180
                  if( fos_angles(ia,0) > 0  &&  fabs(fos_angles(ia,0)) < 180 &&
                      fos_angles(ia,1) != 0  && fabs(fos_angles(ia,1)) < 180 )
                    {
                      // Stretch factor (a max of 100 is applied)
                      const Numeric stretch = max( 100.0, 
                                     1.0/fabs(cos(DEG2RAD*fos_angles(ia,1))) );
                      const Numeric lat0 = rte_pos[1];
                      for( Index i=0; i<lat_grid.nelem(); i++ )
                        { 
                          lat_stretched[i] = lat0 + stretch*(lat_grid[i]-lat0); 
                        }
                    }
                  else
                    { lat_stretched = lat_grid; }
                  
                  iy_clearsky_agendaExecute( ws, tmp, iy_error, iy_error_type,
                                             iy_aux, diy_dx, 0, 
                                             rte_pos, rte_los,
                                             iy_transmission, 0, jacobian_do, 
                                             iy_aux_do, f_grid, p_grid, 
                                             lat_stretched, lon_grid, t_field, 
                                             vmr_field, iy_clearsky_agenda );
                  fos_y(ia,joker,joker) = tmp;
                }
            }
          
        }
      else if( atmosphere_dim == 3 )
        {
          for( Index ia=0; ia<nfosa; ia++ )
            { 
              iy_clearsky_agendaExecute( ws, tmp, iy_error, iy_error_type,
                                         iy_aux, diy_dx, 0, 
                                         rte_pos, fos_angles(ia,Range(0,2)),
                                         iy_transmission, 0, jacobian_do, 
                                         iy_aux_do, f_grid, p_grid, 
                                         lat_grid, lon_grid, t_field, 
                                         vmr_field, iy_clearsky_agenda );
              fos_y(ia,joker,joker) = tmp;
            }
        }
    }
  else
    {
      // The azimuth information is lost for 1D and 2D. Then not possible to
      // handle the latitude stretching for 2D. In fact, the latitude grid
      // should be compressed for many angle combinations!
      if( atmosphere_dim == 2 )
        throw runtime_error( 
             "For atmosphere_dim = 2, only single scattering can be used." );

      Index nlos = 1 + (atmosphere_dim==3);

      for( Index ia=0; ia<nfosa; ia++ )
        { 
          iyFOS( ws, tmp, iy_error, iy_error_type, iy_aux, diy_dx,
                 rte_pos, fos_angles(ia,Range(0,nlos)), iy_aux_do, jacobian_do,
                 atmosphere_dim, p_grid, lat_grid, lon_grid, 
                 z_field, t_field, vmr_field, r_geoid, z_surface, 
                 cloudbox_on, cloudbox_limits, stokes_dim, f_grid, 
                 ppath_step_agenda, emission_agenda, 
                 abs_scalar_gas_agenda, iy_clearsky_agenda, iy_transmission, 
                 pnd_field, scat_data_raw, opt_prop_gas_agenda, fos_y_agenda, 
                 fos_angles, fos_use_mean_scat_data, fos_n, fos_i+1 );

          fos_y(ia,joker,joker) = tmp;
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyFOS(
        Workspace&                     ws,
        Matrix&                        iy,
        Matrix&                        iy_error,
        Index&                         iy_error_type,
        Tensor3&                       iy_aux,
        ArrayOfTensor3&                diy_dx,
  const Vector&                        rte_pos,
  const Vector&                        rte_los,
  const Index&                         iy_aux_do,
  const Index&                         jacobian_do,
  const Index&                         atmosphere_dim,
  const Vector&                        p_grid,
  const Vector&                        lat_grid,
  const Vector&                        lon_grid,
  const Tensor3&                       z_field,
  const Tensor3&                       t_field,
  const Tensor4&                       vmr_field,
  const Matrix&                        r_geoid,
  const Matrix&                        z_surface,
  const Index&                         cloudbox_on,
  const ArrayOfIndex&                  cloudbox_limits,
  const Index&                         stokes_dim,
  const Vector&                        f_grid,
  const Agenda&                        ppath_step_agenda,
  const Agenda&                        emission_agenda,
  const Agenda&                        abs_scalar_gas_agenda,
  const Agenda&                        iy_clearsky_agenda,
  const Tensor3&                       iy_transmission,
  const Tensor4&                       pnd_field,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                        opt_prop_gas_agenda,
  const Agenda&                        fos_y_agenda,
  const Matrix&                        fos_angles,
  const Index&                         fos_use_mean_scat_data,
  const Index&                         fos_n,
  const Index&                         fos_i )
{
  // Input checks
  if( !cloudbox_on )
    throw runtime_error( "The cloudbox must be defined to use this method." );
  if( fos_angles.ncols() != 3 )
    throw runtime_error( "The WSV *fos_angles* must have three columns." );
  if( max(fos_angles) <= PI )
    throw runtime_error( 
                  "The WSV *fos_angles* shall be in degrees (not radians)." );
  if( min(fos_angles(joker,0))<0 || max(fos_angles(joker,0))>180 )
    throw runtime_error( 
                 "The zenith angles in *fos_angles* shall be inside [0,180]." );
  if( min(fos_angles(joker,1))<-180 || max(fos_angles(joker,1))>180 )
    throw runtime_error( 
             "The azimuth angles in *fos_angles* shall be inside [-180,180]." );
  if( fos_angles(joker,2).sum() < 6 || fos_angles(joker,2).sum() > 20 )
    throw runtime_error( 
     "The sum of integration weights in *fos_angles* shall be inside [2,20]." );
  if( fos_n < 0 )
    throw runtime_error( "The WSV *fos_n* must be >= 0." );
  if( fos_i < 0 )
    throw runtime_error( "The WSV *fos_i* must be >= 0." );

  // Determine ppath through the cloudbox
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 0 );

  // Check radiative background
  const Index bkgr = ppath_what_background( ppath );
  if( bkgr == 2 )
    throw runtime_error( "Observations where (unscattered) propagation path "
                         "hits the surface are not yet handled by FOS." );
  assert( bkgr == 3 );

  // Get iy for unscattered direction
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector rte_pos2, rte_los2;
  rte_pos2 = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
  rte_los2 = ppath.los(ppath.np-1,joker);
  //
  iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                             iy_aux, diy_dx, 0, rte_pos2, rte_los2,
                             iy_transmission, 0, jacobian_do, iy_aux_do, 
                             f_grid, p_grid, lat_grid, lon_grid, t_field, 
                             vmr_field, iy_clearsky_agenda );

  // RT for part inside cloudbox
  //
  const Index   np  = ppath.np;
  //
  if( np > 1 )
    {
      // Some sizes
      const Index   nf   = f_grid.nelem();
      const Index   nabs = vmr_field.nbooks();
      const Index   npar = pnd_field.nbooks();

      // Containers for atmospheric data along the LOS
      Vector   ppath_t(np), ppath_p(np);
      Matrix   ppath_vmr(nabs,np);
      Matrix   ppath_pnd(npar,np);

      // Pressure:
      {
        Matrix   itw_p(np,2);
        interpweights( itw_p, ppath.gp_p );
        itw2p( ppath_p, p_grid, ppath.gp_p, itw_p );
      }
      
      // Temperature and VMR
      {
        Matrix   itw_field;
        //
        interp_atmfield_gp2itw( itw_field, atmosphere_dim, 
                                ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
        interp_atmfield_by_itw( ppath_t,  atmosphere_dim, t_field, ppath.gp_p, 
                                ppath.gp_lat, ppath.gp_lon, itw_field );
        for( Index is=0; is<nabs; is++ )
          {
            interp_atmfield_by_itw( ppath_vmr(is,joker), atmosphere_dim, 
                                    vmr_field(is,joker,joker,joker), 
                                    ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                                    itw_field );
          }
      }

      // PND
      {
        Matrix   itw_field;
        //
        interp_cloudfield_gp2itw( itw_field, ppath.gp_p, ppath.gp_lat, 
                               ppath.gp_lon, atmosphere_dim, cloudbox_limits );
        for( Index ip=0; ip<npar; ip++ )
          {
            interp_atmfield_by_itw( ppath_pnd(ip,joker), atmosphere_dim,
                                   pnd_field(ip,joker,joker,joker), ppath.gp_p, 
                                   ppath.gp_lat, ppath.gp_lon, itw_field );
          }
      }

      // Scattering properties
      Array<ArrayOfSingleScatteringData>   scat_data;
      Index   nfs, ivf;
      //
      if( fos_use_mean_scat_data )
        {
          nfs = 1;  ivf = 0;
          scat_data.resize( 1 );
          Vector f_dummy( 1, (f_grid[0]+f_grid[nf-1])/2.0 );
          scat_data_monoCalc( scat_data[0], scat_data_raw, f_dummy, 0 );
        }
      else
        {
          nfs = nf; ivf = 1;
          scat_data.resize( nf );
          for( Index iv=0; iv<nf; iv++ )
            { scat_data_monoCalc( scat_data[iv], scat_data_raw, f_grid, iv ); }
        }

      // Number of FOS angles
      const Index   nfosa = fos_angles.nrows();

      // Data for end points of path step
      // (1 first in propagation direction. Variable 2 for quantities not
      // always calculated must be initiated to 0 (i.e. clear-sky).)
      Matrix asg1(nf,nabs), asg2(nf, nabs);       // Absorption scalar gas
      Vector b1(nf), b2(nf);                      // Emission source function
      bool   cs1, cs2=true;                       // True if no particles 
      Tensor3 emg1(nf,stokes_dim,stokes_dim);     // Extinction matrix
      Tensor3 emg2(nf,stokes_dim,stokes_dim,0.0);
      Matrix  avg1(nf,stokes_dim,0.0);            // Gas absorption vector
      Matrix  avg2(nf,stokes_dim,0.0);
      Tensor3 emp1(nfs,stokes_dim,stokes_dim);    // Particle exinction matrix
      Tensor3 emp2(nfs,stokes_dim,stokes_dim,0.0);
      Matrix  avp1(nfs,stokes_dim);               // Particle absorption vector
      Matrix  avp2(nfs,stokes_dim,0.0); 
      Matrix  s1(nf,stokes_dim);                  // Scattering source term
      Matrix  s2(nf,stokes_dim,0.0);

      // Loop ppath steps
      for( Index ip=np-1; ip>=0; ip-- )
        {
          // Move data one step (new 1 is old 2)
          asg1 = asg2;    
          b1   = b2;
          cs1  = cs2;
          emg1 = emg2;
          avg1 = avg2;
          emp1 = emp2;
          avp1 = avp2;
          s1   = s2;

          // Gas absorption and blackbody radiation
          abs_scalar_gas_agendaExecute( ws, asg2, -1, ppath_p[ip], ppath_t[ip], 
                                  ppath_vmr(joker,ip), abs_scalar_gas_agenda );
          emission_agendaExecute( ws, b2, ppath_t[ip], emission_agenda );

          // Extinction matrix and absorption vector for gases
          // (Not needed for clear-sky, but must anyhow be calculated to get
          // the shifting from point 2 to point 1 correct.
          opt_prop_gas_agendaExecute( ws, emg2, avg2, -1, asg2, 
                                                         opt_prop_gas_agenda );
              
          // Any particles?
          if( max(ppath_pnd(joker,ip)) < 1e-3 ) { cs2 = true; }
          else                                  { cs2 = false; }

          // If clear-sky, set all scattering quantities to zero. 
          // Otherwise make calculations.
          if( cs2 )
            {
              emp2 = 0.0;
              avp2 = 0.0;
              s2   = 0;
            }
          else
            {
              // Direction of outgoing scattered radiation (which is 
              // reversed to LOS).
              // A viewing direction of aa=0 is assumed for 1D. This
              // corresponds to positive za for 2D. Note that rte_los2 is 
              // only used for extracting scattering properties. 
              //
              rte_los2.resize(2);
              //
              Numeric za = ppath.los(ip,0);
              //
              if( atmosphere_dim == 1 )
                { 
                  rte_los2[0] = 180 - za; 
                  rte_los2[1] = 180; 
                }
              else if( atmosphere_dim == 2 )
                {
                  rte_los2[0] = 180 - fabs( za ); 
                  if( za >= 0 )
                    { rte_los2[1] = 180; }
                  else
                    { rte_los2[1] = 0; }
                }
              else if( atmosphere_dim == 3 )
                { 
                  rte_los2[0] = 180 - za; 
                  rte_los2[1] = ppath.los(ip,1) + 180; 
                  if( rte_los2[1] > 180 )
                    { rte_los2[1] -= 360; }
                }

              // Particle properties
              for( Index iv=0; iv<nfs; iv++ )
                { 
                  Matrix dummy1( stokes_dim, stokes_dim );
                  Vector dummy2( stokes_dim );
                  opt_propCalc( dummy1, dummy2, 
                                rte_los2[0], rte_los2[1], scat_data[iv], 
                                stokes_dim, ppath_pnd(joker,ip), ppath_t[ip] );
                  emp2(iv,joker,joker) = dummy1;
                  avp2(iv,joker)       = dummy2;
                }

              // Determine incoming radiation (here Y, WSV is fos_y)
              Tensor3  Y;
              fos_y_agendaExecute( ws, Y, rte_pos2, fos_angles, fos_n, fos_i, 
                                                                fos_y_agenda );

              // Determine phase matrix for fov_angles
              Tensor4  P( nfosa, nfs, stokes_dim, stokes_dim );
              Matrix   P1( stokes_dim, stokes_dim );
              //
              for( Index ia=0; ia<nfosa; ia++ )
                { 
                  for( Index iv=0; iv<nfs; iv++ )
                    {
                      pha_mat_singleCalc( P1, rte_los2[0], rte_los2[1],
                                          fos_angles(ia,0), fos_angles(ia,1),
                                          scat_data[iv], stokes_dim, 
                                          ppath_pnd(joker,ip), ppath_t[ip] );
                      P(ia,iv,joker,joker) = P1;
                    }
                }

              // Scattering source term
              s2 = 0.0;
              for( Index iv=0; iv<nf; iv++ )
                { 
                  Vector sp(stokes_dim);
                  for( Index ia=0; ia<nfosa; ia++ )
                    { 
                      mult( sp, P(ia,iv*ivf,joker,joker), Y(ia,iv,joker) );
                      sp           *= fos_angles(ia,2);
                      s2(iv,joker) += sp;
                    }
                }
            }

          // RT of ppath step (nothing to do when at upper point)
          if( ip < np-1 )
            {
              // Clearsky step
              if( ( cs1 && cs2 )  ||  fos_n == 0 )
                {
                  // Loop frequencies
                  for( Index iv=0; iv<nf; iv++ )
                    { 
                      const Numeric step_tr = exp( -ppath.l_step[ip] * 0.5 *
                             ( asg1(iv,joker).sum() + asg2(iv,joker).sum() ) );
                      //
                      iy(iv,0)  = iy(iv,0) * step_tr + (1-step_tr) * 0.5 *
                                                  (b1[iv] + b2[iv]);
                      //
                      for( Index is=1; is<stokes_dim; is++ )
                        { iy(iv,is) *= step_tr; }
                    }
                }
          
              // RT step with scattering
              else
                {
                  // Loop frequencies
                  for( Index iv=0; iv<nf; iv++ )
                    {
                      // Calculate average of absorption, extinction etc.
                      Matrix ext_mat(stokes_dim,stokes_dim,0.0);
                      Vector abs_vec(stokes_dim,0.0), s(stokes_dim,0.0);
                      Numeric b = 0.5 * ( b1[iv] + b2[iv] );
                      //
                      for( Index i1=0; i1<stokes_dim; i1++ )
                        {
                          abs_vec[i1] = 0.5 * ( avg1(iv,i1) + avg2(iv,i1) +
                                        avp1(iv*ivf,i1) + avp2(iv*ivf,i1) ); 
                          s[i1]  = 0.5 * ( s1(iv,i1) + s2(iv,i1) );
                          //
                          for( Index i2=0; i2<stokes_dim; i2++ )
                            {
                              ext_mat(i1,i2) = 0.5 * ( 
                                     emg1(iv,    i1,i2) + emg2(iv,    i1,i2) + 
                                     emp1(iv*ivf,i1,i2) + emp2(iv*ivf,i1,i2) );
                            }
                        }

                      // RT for step
                      Matrix trans_mat(stokes_dim,stokes_dim); 
                      //
                      rte_step_std( iy(iv,joker), trans_mat, ext_mat, abs_vec,
                                                     s, ppath.l_step[ip], b );
                    }
                }
            }
        } 
    }
}


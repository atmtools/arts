/* Copyright (C) 2012
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
void fos_yStandard(Workspace&          ws,
                   Tensor3&            fos_y,
                   Matrix&             iy_error,
                   Index&              iy_error_type,
                   Matrix&             iy_aux,
                   ArrayOfTensor3&     diy_dx,
                   const Vector&       rte_pos,
                   const Index&        atmosphere_dim,
                   const Vector&       p_grid,
                   const Vector&       lat_grid,
                   const Tensor3&      z_field,
                   const Tensor3&      t_field,
                   const Tensor4&      vmr_field,
                   const Tensor3&      edensity_field,
                   const Index&        cloudbox_on,
                   const ArrayOfIndex& cloudbox_limits,
                   const Index&        stokes_dim,
                   const Vector&       f_grid,
                   const Agenda&       ppath_agenda,
                   const Agenda&       emission_agenda,
                   const Agenda&       abs_scalar_gas_agenda,
                   const Agenda&       iy_clearsky_agenda,
                   const Tensor3&      iy_transmission,
                   const Tensor4&      pnd_field,
                   const ArrayOfSingleScatteringData&   scat_data_raw,
                   const Agenda&       opt_prop_gas_agenda,
                   const Agenda&       fos_y_agenda,
                   const Matrix&       fos_angles,
                   const Index&        use_mean_scat_data,
                   const Index&        fos_n,
                   const Index&        fos_i,
                   const Verbosity&    verbosity)
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
                                             iy_aux, diy_dx, 0, iy_transmission,
                                             rte_pos, fos_angles(ia,Range(0,1)),
                                             0, jacobian_do, t_field, z_field, 
                                             vmr_field, -1,
                                             iy_clearsky_agenda );
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
                  
                  throw runtime_error( "2D FOS can not be used presently "
                       "as the \"lat-stretch\" approach not can be handled!" );

                  iy_clearsky_agendaExecute( ws, tmp, iy_error, iy_error_type,
                                             iy_aux, diy_dx, 0, iy_transmission,
                                             rte_pos, rte_los, 0, jacobian_do, 
                                             t_field, z_field, vmr_field, -1, 
                                             iy_clearsky_agenda );
                  fos_y(ia,joker,joker) = tmp;
                }
            }
          
        }
      else if( atmosphere_dim == 3 )
        {
          for( Index ia=0; ia<nfosa; ia++ )
            { 
              iy_clearsky_agendaExecute( ws, tmp, iy_error, iy_error_type,
                                         iy_aux, diy_dx, 0, iy_transmission, 
                                         rte_pos, fos_angles(ia,Range(0,2)),
                                         0, jacobian_do, t_field, z_field, 
                                         vmr_field, -1, iy_clearsky_agenda );
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
                 iy_transmission, rte_pos, fos_angles(ia,Range(0,nlos)), 
                 jacobian_do, atmosphere_dim, p_grid, 
                 z_field, t_field, vmr_field, edensity_field,
                 cloudbox_on, cloudbox_limits, stokes_dim, f_grid, 
                 ppath_agenda, emission_agenda, 
                 abs_scalar_gas_agenda, iy_clearsky_agenda, 
                 pnd_field, scat_data_raw, opt_prop_gas_agenda, fos_y_agenda, 
                 fos_angles, use_mean_scat_data, fos_n, fos_i+1, verbosity);

          fos_y(ia,joker,joker) = tmp;
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyFOS(Workspace&          ws,
           Matrix&             iy,
           Matrix&             iy_error,
           Index&              iy_error_type,
           Matrix&             iy_aux,
           ArrayOfTensor3&     diy_dx,
           const Tensor3&      iy_transmission,
           const Vector&       rte_pos,
           const Vector&       rte_los,
           const Index&        jacobian_do,
           const Index&        atmosphere_dim,
           const Vector&       p_grid,
           const Tensor3&      z_field,
           const Tensor3&      t_field,
           const Tensor4&      vmr_field,
           const Tensor3&      edensity_field,
           const Index&        cloudbox_on,
           const ArrayOfIndex& cloudbox_limits,
           const Index&        stokes_dim,
           const Vector&       f_grid,
           const Agenda&       ppath_agenda,
           const Agenda&       emission_agenda,
           const Agenda&       abs_scalar_gas_agenda,
           const Agenda&       iy_clearsky_agenda,
           const Tensor4&      pnd_field,
           const ArrayOfSingleScatteringData&   scat_data_raw,
           const Agenda&       opt_prop_gas_agenda,
           const Agenda&       fos_y_agenda,
           const Matrix&       fos_angles,
           const Index&        use_mean_scat_data,
           const Index&        fos_n,
           const Index&        fos_i,
           const Verbosity&    verbosity)
{
  // Input checks
  if( jacobian_do )
    throw runtime_error( 
     "This method does not yet provide any jacobians (jacobian_do must be 0)" );
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
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 1, -1,
                       t_field, z_field, vmr_field, edensity_field, -1,
                       ppath_agenda );

  // Check radiative background
  const Index bkgr = ppath_what_background( ppath );
  if( bkgr == 2 )
    throw runtime_error( "Observations where (unscattered) propagation path "
                         "hits the surface inside the cloudbox are not yet "
                         "handled by this method." );
  assert( bkgr == 3 );

  // Get atmospheric and RT quantities for each ppath point/step (inside box)
  // 
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
  Matrix    ppath_vmr, ppath_pnd;
  // "rtvars"
  Matrix    ppath_emission, ppath_tau;
  Tensor3   wind_field_dummy(0,0,0), iy_trans_new;
  Tensor3   ppath_asp_abs_vec, ppath_pnd_abs_vec, total_transmission;
  Tensor4   ppath_asp_ext_mat, ppath_pnd_ext_mat, ppath_transmission;
  Array<ArrayOfSingleScatteringData>  scat_data;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_field_dummy, wind_field_dummy, wind_field_dummy );

      // Particle number densities
      get_ppath_pnd( ppath_pnd, 
                     ppath, atmosphere_dim, cloudbox_limits, pnd_field );

      // Absorption and optical thickness for each step
      get_ppath_cloudrtvars( ws, ppath_asp_abs_vec, ppath_asp_ext_mat,
                            ppath_pnd_abs_vec, ppath_pnd_ext_mat, ppath_transmission, 
                            total_transmission, ppath_emission, scat_data,
                            abs_scalar_gas_agenda, emission_agenda, opt_prop_gas_agenda,
                            ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                            ppath_wind_v, ppath_wind_w, ppath_pnd, use_mean_scat_data,
                            scat_data_raw, stokes_dim, f_grid, atmosphere_dim, 1,
                            verbosity);
    }
  else // Just in case, should not happen
    { assert( 0 ); }

  // iy_transmission
  //
  iy_transmission_mult( iy_trans_new, iy_transmission, total_transmission );

  // Get iy for unscattered direction
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector  rte_pos2;
  {
    Vector  rte_los2;
    rte_pos2 = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
    rte_los2 = ppath.los(ppath.np-1,joker);
    //
    iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                               iy_aux, diy_dx, 0, iy_trans_new,
                               rte_pos2, rte_los2, 0, jacobian_do, t_field, 
                               z_field, vmr_field, -1, iy_clearsky_agenda );
  }

  // RT for part inside cloudbox
  //
  if( np > 1 )
    {
      // General variables
      const Index   nf    = f_grid.nelem();
      const Index   nfosa = fos_angles.nrows();   // Number of FOS angles
      Matrix  s1(nf,stokes_dim);                  // Scattering source term
      Matrix  s2(nf,stokes_dim,0.0);

      // Help variables for handling of *use_mean_scat_data*
      Index   nfs, ivf;
      //
      if( use_mean_scat_data )
        { nfs = 1;  ivf = 0; }
      else
        { nfs = nf; ivf = 1; }

      // Loop ppath steps
      for( Index ip=np-1; ip>=0; ip-- )
        {              
          // Update scattering source term (new 1 is old 2)
          s1 = s2;

          // Scattering source term (is zero if no particles)
          if( max(ppath_pnd(joker,ip)) < 1e-3 )
            {
              s2 = 0.0;
            }
          else
            {
              // Determine incoming radiation (here Y, WSV equals fos_y)
              Tensor3  Y;
              fos_y_agendaExecute( ws, Y, rte_pos2, fos_angles, fos_n, fos_i, 
                                                                fos_y_agenda );


              // Direction of outgoing scattered radiation (which is reversed
              // to LOS). Note that this rte_los2 is only used for extracting
              // scattering properties.
              Vector rte_los2;
              mirror_los( rte_los2, ppath.los(ip,joker), atmosphere_dim );

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
                                          ppath_pnd(joker,ip), ppath_t[ip],
                                          verbosity );
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
              // Loop frequencies
              for( Index iv=0; iv<nf; iv++ )
                {
                  // Calculate average of absorption, extinction etc.
                  Matrix  ext_mat( stokes_dim, stokes_dim,0 );
                  Vector  abs_vec(stokes_dim,0.0);
                  Vector  s(stokes_dim,0.0);
                  Numeric b = 0.5 * ( ppath_emission(iv,ip) + 
                                      ppath_emission(iv,ip+1) );
                  for( Index is1=0; is1<stokes_dim; is1++ )
                    { 
                      s[is1]  = 0.5 * ( s1(iv,is1) + s2(iv,is1) );
                      abs_vec[is1] = 0.5 * ( 
                                          ppath_asp_abs_vec(iv,is1,ip+1) +
                                          ppath_asp_abs_vec(iv,is1,ip)   +
                                          ppath_pnd_abs_vec(iv,is1,ip+1) +
                                          ppath_pnd_abs_vec(iv,is1,ip) );
                      for( Index is2=0; is2<stokes_dim; is2++ )
                        {
                          ext_mat(is1,is2) = 0.5 * (
                                          ppath_asp_ext_mat(iv,is1,is2,ip+1) +
                                          ppath_asp_ext_mat(iv,is1,is2,ip)   +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip+1) +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip) );
                        }
                    }

                  // RT for step
                  Matrix trans_mat(stokes_dim,stokes_dim); 
                  //
                  rte_step_std( iy(iv,joker), trans_mat, ext_mat, abs_vec,
                                                       s, ppath.lstep[ip], b );
                }
            }
        } 
    }
}


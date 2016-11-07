/* copyright (C) 2003-2012 Cory Davis <cory.davis@metservice.com>
                            
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
  \file   montecarlo.cc
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2003-06-19 

  \brief  functions used by MCGeneral
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <sstream>

#include "auto_md.h"
#include "geodetic.h"
#include "montecarlo.h"
#include "mc_interp.h"

//! clear_rt_vars_at_gp
/*! 

  Calculates a bunch of atmospheric variables at the end of a ppath.
   
\author Cory Davis
\date 2005-02-19?
*/  
void clear_rt_vars_at_gp(Workspace&          ws,
                         MatrixView          ext_mat_mono,
                         VectorView          abs_vec_mono,
                         Numeric&            temperature,
                         const Agenda&       propmat_clearsky_agenda,
                         const Numeric&      f_mono,
                         const GridPos&      gp_p,
                         const GridPos&      gp_lat,
                         const GridPos&      gp_lon,
                         ConstVectorView     p_grid,
                         ConstTensor3View    t_field,
                         ConstTensor4View    vmr_field)
{
  const Index   ns = vmr_field.nbooks();
  Vector t_vec(1);   //vectors are required by interp_atmfield_gp2itw etc.
  Vector   p_vec(1); //may not be efficient with unecessary vectors
  Matrix itw_p(1,2);
  ArrayOfGridPos ao_gp_p(1),ao_gp_lat(1),ao_gp_lon(1);
  Matrix   vmr_mat(ns,1), itw_field;
  
  //local versions of workspace variables
  Matrix  local_abs_vec;
  Tensor3 local_ext_mat, local_nlte_source_dummy;
  Tensor4 local_propmat_clearsky;
  ArrayOfTensor3 local_partial_dummy; // This is right since there should be only clearsky partials
  ArrayOfMatrix local_dnlte_dx_source_dummy,local_nlte_dsource_dx_dummy;
  ao_gp_p[0]=gp_p;
  ao_gp_lat[0]=gp_lat;
  ao_gp_lon[0]=gp_lon;
  
  // Determine the pressure 
  interpweights( itw_p, ao_gp_p );
  itw2p( p_vec, p_grid, ao_gp_p, itw_p );
  
  // Determine the atmospheric temperature and species VMR 
  //
  interp_atmfield_gp2itw( itw_field, 3, ao_gp_p, ao_gp_lat, ao_gp_lon );
  //
  interp_atmfield_by_itw( t_vec,  3, t_field, ao_gp_p, ao_gp_lat, ao_gp_lon, 
                          itw_field );
  // 
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( vmr_mat(is, joker), 3, 
                              vmr_field(is, joker, joker, joker), 
                              ao_gp_p, ao_gp_lat, 
                              ao_gp_lon, itw_field );
    }


  temperature = t_vec[0];

  const Vector rtp_mag_dummy(3,0);
  const Vector ppath_los_dummy;
  const Vector temperature_nlte_dummy(0);
  
  //calcualte absorption coefficient
  propmat_clearsky_agendaExecute(ws, local_propmat_clearsky,local_nlte_source_dummy,local_partial_dummy,local_dnlte_dx_source_dummy,local_nlte_dsource_dx_dummy,
                                 ArrayOfRetrievalQuantity(0), Vector(1, f_mono), rtp_mag_dummy,
                                 ppath_los_dummy,p_vec[0],
                                 temperature, temperature_nlte_dummy,
                                 vmr_mat(joker, 0),
                                 propmat_clearsky_agenda);

  opt_prop_sum_propmat_clearsky(local_ext_mat, local_abs_vec, local_propmat_clearsky);
  
  ext_mat_mono=local_ext_mat(0, joker, joker);
  abs_vec_mono=local_abs_vec(0,joker);
}



//! cloudy_rt_vars_at_gp
/*! 

  Calculates a bunch of atmospheric variables at the end of a ppath.
   
\author Cory Davis
\date 2005-02-19?
*/  


void cloudy_rt_vars_at_gp(Workspace&           ws,
                          MatrixView           ext_mat_mono,
                          VectorView           abs_vec_mono,
                          VectorView           pnd_vec,
                          Numeric&             temperature,
                          const Agenda&        propmat_clearsky_agenda,
                          const Index          stokes_dim,
                          const Numeric&       f_mono,
                          const GridPos&       gp_p,
                          const GridPos&       gp_lat,
                          const GridPos&       gp_lon,
                          ConstVectorView      p_grid_cloud,
                          ConstTensor3View     t_field_cloud,
                          ConstTensor4View     vmr_field_cloud,
                          const Tensor4&       pnd_field,
                          const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                          const ArrayOfIndex&  cloudbox_limits,
                          const Vector&        rte_los,
                          const Verbosity&     verbosity
                          )

{
  const Index   ns = vmr_field_cloud.nbooks();
  const Index N_se = pnd_field.nbooks();
  Matrix  pnd_ppath(N_se,1);
  Vector t_ppath(1);
  Vector   p_ppath(1);//may not be efficient with unecessary vectors
  Vector temperature_nlte_dummy(0);
  Matrix itw_p(1,2);
  ArrayOfGridPos ao_gp_p(1),ao_gp_lat(1),ao_gp_lon(1);
  Matrix   vmr_ppath(ns,1), itw_field;
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);
  Vector abs_vec_part(stokes_dim, 0.0);

  //local versions of workspace variables
  ArrayOfTensor3 local_partial_dummy; // This is right since there should be only clearsky partials
  ArrayOfMatrix local_dnlte_dx_source_dummy; // This is right since there should be only clearsky partials
  ArrayOfMatrix local_nlte_dsource_dx_dummy; // This is right since there should be only clearsky partials
  Tensor4 local_propmat_clearsky;
  Tensor3 local_nlte_source_dummy;//FIXME: Do this right?
  Matrix  local_abs_vec;
  Tensor3 local_ext_mat;

  ao_gp_p[0]=gp_p;
  ao_gp_lat[0]=gp_lat;
  ao_gp_lon[0]=gp_lon;
  

  cloud_atm_vars_by_gp(p_ppath,t_ppath,vmr_ppath,pnd_ppath,ao_gp_p,
                       ao_gp_lat,ao_gp_lon,cloudbox_limits,p_grid_cloud,
                       t_field_cloud, vmr_field_cloud,pnd_field);
   
  temperature = t_ppath[0];

  const Vector rtp_mag_dummy(3,0);
  const Vector ppath_los_dummy;
  
  //rtp_vmr    = vmr_ppath(joker,0);
  propmat_clearsky_agendaExecute(ws, local_propmat_clearsky,
                                 local_nlte_source_dummy,local_partial_dummy,
                                 local_dnlte_dx_source_dummy,local_nlte_dsource_dx_dummy,
                                 ArrayOfRetrievalQuantity(0),
                                 Vector(1, f_mono), rtp_mag_dummy,
                                 ppath_los_dummy,p_ppath[0],
                                 temperature, temperature_nlte_dummy,
                                 vmr_ppath(joker, 0),
                                 propmat_clearsky_agenda);
  
  opt_prop_sum_propmat_clearsky(local_ext_mat, local_abs_vec, local_propmat_clearsky);

  
  ext_mat_mono=local_ext_mat(0, joker, joker);
  abs_vec_mono=local_abs_vec(0,joker);
  ext_mat_part=0.0;
  abs_vec_part=0.0;

  Vector sca_dir;
  mirror_los( sca_dir, rte_los, 3 );
  //
  pnd_vec = pnd_ppath(joker, 0);
  opt_propCalc(ext_mat_part,abs_vec_part,sca_dir[0],sca_dir[1],scat_data_mono,
               stokes_dim, pnd_vec, temperature,verbosity);
  
  ext_mat_mono += ext_mat_part;
  abs_vec_mono += abs_vec_part;

}



//! cloud_atm_vars_by_gp
/*! 

  Returns pressure, temperature, VMRs and PNDs, at points corresponding
  to arrays of gridpositions gp_p, gp_lat, and gp_lon.  The field and grid 
  input variables all span only the cloudbox

  \param pressure  Output: a vector of pressures
  \param temperature  Output: a vector of temperatures
  \param vmr          Output: a n_species by n_p matrix of VMRs
  \param pnd          Output: a n_scatelem by n_p matrix of PNDs
  \param gp_p         an array of pressure gridpoints
  \param gp_lat       an array of latitude gridpoints
  \param gp_lon       an array of longitude gridpoints
  \param cloudbox_limits  the WSV
  \param p_grid_cloud the subset of the p_grid corresponding to the cloudbox
  \param t_field_cloud  the t_field within the cloudbox
  \param vmr_field_cloud the t_field within the cloudbox
  \param pnd_field             The WSV
\author Cory Davis
\date 2005-06-07
*/

void cloud_atm_vars_by_gp(
                          VectorView pressure,
                          VectorView temperature,
                          MatrixView vmr,
                          MatrixView pnd,
                          const ArrayOfGridPos& gp_p,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const ArrayOfIndex&   cloudbox_limits,
                          ConstVectorView    p_grid_cloud,
                          ConstTensor3View   t_field_cloud,
                          ConstTensor4View   vmr_field_cloud,
                          ConstTensor4View   pnd_field
                          )
{
  Index np=gp_p.nelem();
  assert(pressure.nelem()==np);
  Index ns=vmr_field_cloud.nbooks();
  Index N_se=pnd_field.nbooks();
  ArrayOfGridPos gp_p_cloud   = gp_p;
  ArrayOfGridPos gp_lat_cloud = gp_lat;
  ArrayOfGridPos gp_lon_cloud = gp_lon;
  Index atmosphere_dim=3;

  for (Index i = 0; i < np; i++ ) 
    {
      // Calculate grid positions inside the cloud. 
      gp_p_cloud[i].idx -= cloudbox_limits[0];
      gp_lat_cloud[i].idx -= cloudbox_limits[2];
      gp_lon_cloud[i].idx -= cloudbox_limits[4];
    }      

  const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
  const Index n2 = cloudbox_limits[3] - cloudbox_limits[2];
  const Index n3 = cloudbox_limits[5] - cloudbox_limits[4];
  gridpos_upperend_check( gp_p_cloud[0],      n1 );
  gridpos_upperend_check( gp_p_cloud[np-1],   n1 );
  gridpos_upperend_check( gp_lat_cloud[0],    n2 );
  gridpos_upperend_check( gp_lat_cloud[np-1], n2 );
  gridpos_upperend_check( gp_lon_cloud[0],    n3 );
  gridpos_upperend_check( gp_lon_cloud[np-1], n3 );
  
  // Determine the pressure at each propagation path point
  Matrix   itw_p(np,2);
  //
  //interpweights( itw_p, ppath.gp_p );      
  interpweights( itw_p, gp_p_cloud );
  itw2p( pressure, p_grid_cloud, gp_p_cloud, itw_p );
  
  // Determine the atmospheric temperature and species VMR at 
  // each propagation path point
  Matrix   itw_field;
  //
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, 
                          gp_p_cloud, gp_lat_cloud, gp_lon_cloud );
  //
  interp_atmfield_by_itw( temperature,  atmosphere_dim, t_field_cloud, 
                          gp_p_cloud, gp_lat_cloud, gp_lon_cloud, 
                          itw_field );
  // 
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( vmr(is, joker), atmosphere_dim,  
                              vmr_field_cloud(is, joker, joker, joker), 
                              gp_p_cloud, gp_lat_cloud, 
                              gp_lon_cloud, itw_field );
    }
  
  //Determine the particle number density for every scattering element at 
  // each propagation path point
  for( Index i_se=0; i_se<N_se; i_se++ )
    {
      // if grid positions still outside the range the propagation path step 
      // must be outside the cloudbox and pnd is set to zero
      interp_atmfield_by_itw( pnd(i_se, joker), atmosphere_dim,
                              pnd_field(i_se, joker, joker,  joker), 
                              gp_p_cloud, gp_lat_cloud, 
                              gp_lon_cloud, itw_field );
    }
}



//! findZ11max
/*! 
  The direction sampling method requires a bounding value for Z11.
  This returns a vector with the maximum value of Z11 for each scattering
  element.

  \param[out] Z11maxvector Maximum value of Z11 for each scattering element
  \param[in]  scat_data_mono

  \author Cory Davis
  \date   2004-31-1
*/
void findZ11max(Vector& Z11maxvector,
                const ArrayOfArrayOfSingleScatteringData& scat_data_mono)
{
  Z11maxvector.resize(TotalNumberOfElements(scat_data_mono));

  Index i_total = -1;
  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss<scat_data_mono.nelem(); i_ss++)
  {
      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
      {
          i_total++;
          switch(scat_data_mono[i_ss][i_se].ptype){
              case PTYPE_MACROS_ISO:
              {
                  Z11maxvector[i_total]=max(scat_data_mono[i_ss][i_se].pha_mat_data(0,joker,joker,0,0,0,0));
              }
              case PTYPE_HORIZ_AL:
              {
                  Z11maxvector[i_total]=max(scat_data_mono[i_ss][i_se].pha_mat_data(0,joker,joker,0,joker,0,0));
              }
              default:
                  Z11maxvector[i_total]=max(scat_data_mono[i_ss][i_se].pha_mat_data(0,joker,joker,joker,joker,joker,0));
          }
      }
  }
}




//! is_anyptype30
/*!
Some operations in Monte Carlo simulations are different depending on the 
ptype of the scattering elements. This function searches scat_data_mono
to determine if any of the scattering elements have ptype=30.

\author Cory Davis
\date 2004-1-31

*/
bool is_anyptype30(const ArrayOfArrayOfSingleScatteringData& scat_data_mono)
{
    bool anyptype30=false;
    for (Index i_ss = 0; anyptype30 == false && i_ss<scat_data_mono.nelem(); i_ss++)
    {
        for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
        {
            if(scat_data_mono[i_ss][i_se].ptype==PTYPE_HORIZ_AL)
            {
                anyptype30=true;
            }
        }
    }

    return anyptype30;
}



//! mcPathTraceGeneral
/*!
    Performs the tasks of pathlength sampling.

    Ray tracing done (but now only as far as determined by pathlength 
    sampling) and calculation of the evolution operator and several 
    atmospheric variables at the new point.

    The end point of the ray tracing is returned by ppath_step, where the 
    point of concern has index ppath_step.np-1. However, a somehwat dirty trick
    is used here to avoid copying of data. Only ppath.np is adjusted, and
    ppath_step can contain additional points (that should not be used).

    2012-11-15  Patrick Eriksson
    Revised.  Added handling of ppath_step_agenda. Correct handling of ppath
    steps having more than two points.

    2016-10-12  Patrick Eriksson
    Added taustep_limit feature.

    \author Cory Davis
    \date 2005-2-21
*/
void mcPathTraceGeneral(
         Workspace&      ws,
         MatrixView      evol_op,
         Vector&         abs_vec_mono,
         Numeric&        temperature,
         MatrixView      ext_mat_mono,
         Rng&            rng,
         Vector&         rte_pos,
         Vector&         rte_los,
         Vector&         pnd_vec,
         Numeric&        g,
         Ppath&          ppath_step,
         Index&          termination_flag,
         bool&           inside_cloud,
   const Agenda&         ppath_step_agenda,
   const Numeric&        ppath_lmax,
   const Numeric&        ppath_lraytrace,
   const Numeric&        taustep_limit,
   const Agenda&         propmat_clearsky_agenda,
   const Index           stokes_dim,
   const Numeric&        f_mono,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        z_field,
   const Vector&         refellipsoid,
   const Matrix&         z_surface,
   const Tensor3&        t_field,
   const Tensor4&        vmr_field,
   const ArrayOfIndex&   cloudbox_limits,
   const Tensor4&        pnd_field,
   const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
   const Verbosity&      verbosity )
{ 
  ArrayOfMatrix evol_opArray(2); 
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2);
  ArrayOfVector pnd_vecArray(2);
  Matrix        ext_mat(stokes_dim,stokes_dim);
  Matrix        incT(stokes_dim,stokes_dim,0.0);
  Vector        tArray(2);
  Vector        f_grid(1,f_mono);     // Vector version of f_mono
  Matrix        T(stokes_dim,stokes_dim);
  Numeric       k;
  Numeric       ds, dl=-999;
  Index         istep = 0;  // Counter for number of steps
  Matrix        old_evol_op(stokes_dim,stokes_dim);

  CREATE_OUT0;

  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1]=evol_op;

  // range defining cloudbox
  Range p_range(   cloudbox_limits[0], 
                   cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range( cloudbox_limits[2], 
                   cloudbox_limits[3]-cloudbox_limits[2]+1 );
  Range lon_range( cloudbox_limits[4], 
                   cloudbox_limits[5]-cloudbox_limits[4]+1 );

  //initialise Ppath with ppath_start_stepping
  ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, refellipsoid, z_surface,
                        0, cloudbox_limits, false, 
                        rte_pos, rte_los, verbosity );

  // Check if we have already has radiative background
  if( ppath_what_background(ppath_step) )
    { 
      termination_flag = ppath_what_background(ppath_step);
      g = 1;
      return;
    }
  
  // Index in ppath_step of end point considered presently
  Index ip = 0;

  // Is point ip inside the cloudbox?
  inside_cloud = is_gp_inside_cloudbox( ppath_step.gp_p[ip], 
                                        ppath_step.gp_lat[ip], 
                                        ppath_step.gp_lon[ip], 
                                        cloudbox_limits, 0, 3 );

  // Determine radiative properties at point
  if( inside_cloud )
    {
      cloudy_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, pnd_vec, 
                            temperature, propmat_clearsky_agenda, 
                            stokes_dim, f_mono, ppath_step.gp_p[0], 
                            ppath_step.gp_lat[0], ppath_step.gp_lon[0],
                            p_grid[p_range], 
                            t_field(p_range,lat_range,lon_range), 
                            vmr_field(joker,p_range,lat_range,lon_range),
                            pnd_field, scat_data_mono, cloudbox_limits,
                            ppath_step.los(0,joker), verbosity );
    }
  else
    {
      clear_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, temperature, 
                           propmat_clearsky_agenda, f_mono,
                           ppath_step.gp_p[0], ppath_step.gp_lat[0], 
                           ppath_step.gp_lon[0], p_grid, t_field, vmr_field );
      pnd_vec = 0.0;
    }

  // Move the data to end point containers 
  ext_matArray[1] = ext_mat_mono;
  abs_vecArray[1] = abs_vec_mono;
  tArray[1]       = temperature;
  pnd_vecArray[1] = pnd_vec;

  //draw random number to determine end point
  Numeric r = rng.draw();

  termination_flag=0;

  while( (evol_op(0,0)>r) && (!termination_flag) )
    {
      istep++;

      if( istep > 100000 )
        {
          throw runtime_error( "100000 path points have been reached. "
                               "Is this an infinite loop?" );
        }

      evol_opArray[0] = evol_opArray[1];
      ext_matArray[0] = ext_matArray[1];
      abs_vecArray[0] = abs_vecArray[1];
      tArray[0]       = tArray[1];
      pnd_vecArray[0] = pnd_vecArray[1];

      // The algorith to meet taustep_lim:
      // When first calculating a new ppath_step, it assumed that the present
      // ppath position holds the highest extinction. If the extinction at the
      // next position is higher, the criterion is checked and a new ppath_step
      // calculation is triggered if found necessary.
      // This should work in most cases, but is not 100% safe. Consider a case
      // with ppath_lmax = -1 and the extinction is zero at all grid box
      // corners except one. The two "test points" can then both get an
      // extinction of zero, while in fact is non-zero through the grid box and
      // the optical depth is underestimated. But this was handled equally bad
      // before taustep_limit was introduced (2016-10-10, PE) 
      bool  oktaustep = false;
      Index ppath_try = 1;
      const Index lmax_limit = 10;
      
      while( !oktaustep )
        {
          // Shall new ppath_step be calculated?
          if( ip == ppath_step.np-1 ) 
            {
              Numeric lmax = taustep_limit/ext_mat_mono(0,0);
              if( ppath_lmax > 0 )
                { lmax = min( ppath_lmax, lmax ); }
              if( lmax < lmax_limit ) { lmax = lmax_limit; }
              //cout << ppath_try << ", lmax = " << lmax << endl;              
              //Print( ppath_step, 0, verbosity );

              ppath_step_agendaExecute( ws, ppath_step, lmax, ppath_lraytrace,
                                        t_field, z_field, vmr_field, f_grid, 
                                        ppath_step_agenda );
              ip = 1;

              inside_cloud = is_gp_inside_cloudbox( ppath_step.gp_p[ip], 
                                                    ppath_step.gp_lat[ip], 
                                                    ppath_step.gp_lon[ip], 
                                                    cloudbox_limits, 0, 3 );
            }
          else
            { ip++; }

          if( inside_cloud )
            {
              cloudy_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, pnd_vec,
                                    temperature, propmat_clearsky_agenda, 
                                    stokes_dim, f_mono, ppath_step.gp_p[ip],
                                    ppath_step.gp_lat[ip], ppath_step.gp_lon[ip],
                                    p_grid[p_range], 
                                    t_field(p_range,lat_range,lon_range), 
                                    vmr_field(joker,p_range,lat_range,lon_range),
                                    pnd_field, scat_data_mono, cloudbox_limits,
                                    ppath_step.los(ip,joker), verbosity );
            }
          else
            {
              clear_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, temperature, 
                                   propmat_clearsky_agenda, f_mono,
                                   ppath_step.gp_p[ip], ppath_step.gp_lat[ip],
                                   ppath_step.gp_lon[ip], p_grid, t_field, 
                                   vmr_field );
              pnd_vec = 0.0;
            }

          dl = ppath_step.lstep[ip-1];
          
          // Check if taustep_limit criterion is met. OK if:
          // 1. Ppath step already recalculated
          // 2. New ext_mat <= old one
          // 3. New ext_mat bigger, but tau of step still below limit
          if( ppath_try > 1  ||  ext_mat_mono(0,0) <= ext_matArray[0](0,0)  ||  
              (ext_mat_mono(0,0)+ext_matArray[0](0,0))*dl/2 <= taustep_limit )
            {
              oktaustep = true;
            }
          else
            {
              // We trigger a recalculation of ppath_step, from the previous
              // point
              ppath_step.np = ip; 
              ip--;
              ppath_try     = 2;
              // If a background found in first try this has to be reset:
              ppath_set_background( ppath_step, 0 );
            }
        } // while !oktuastep
      
      ext_matArray[1] = ext_mat_mono;
      abs_vecArray[1] = abs_vec_mono;
      tArray[1]       = temperature;
      pnd_vecArray[1] = pnd_vec;
      ext_mat         = ext_matArray[1];
      ext_mat        += ext_matArray[0];  // Factor 2 fixed by using dl/2
      //
      {
        Index extmat_case = 0;
        ext2trans( incT, extmat_case, ext_mat, dl/2 );
      }
      //
      mult( evol_op, evol_opArray[0], incT );
      evol_opArray[1] = evol_op;
     
      if( evol_op(0,0) > r )
        {
          // Check whether hit ground or space.
          // path_step_agenda just detects surface intersections, and
          // if TOA is reached requires a special check.
          // But we are already ready if evol_op<=r
          if( ip == ppath_step.np - 1 )
            {
              if( ppath_what_background(ppath_step) )
                { termination_flag = 2; }   //we have hit the surface
              else if( fractional_gp(ppath_step.gp_p[ip]) >= 
                                          (Numeric)(p_grid.nelem() - 1)-1e-3 )
                { termination_flag = 1; }  //we are at TOA
            }
        }
    } // while


  if( termination_flag ) 
    { //we must have reached the cloudbox boundary
      rte_pos = ppath_step.pos(ip,joker);
      rte_los = ppath_step.los(ip,joker);
      g       = evol_op(0,0);
    }
  else
    {
      //find position...and evol_op..and everything else required at the new
      //scattering/emission point
      // GH 2011-09-14: 
      //   log(incT(0,0)) = log(exp(opt_depth_mat(0, 0))) = opt_depth_mat(0, 0)
      //   Avoid loss of precision, use opt_depth_mat directly
      //k=-log(incT(0,0))/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
      // PE 2013-05-17, Now the above comes directly from ext_mat:
      k  = ext_mat(0,0) / 2;   // Factor 2 as sum of end point values
      ds = log( evol_opArray[0](0,0) / r ) / k;
      g  = k*r;
      Vector x(2,0.0);
      //interpolate atmospheric variables required later
      ArrayOfGridPos gp(1);
      x[1] = dl;
      Vector itw(2);
  
      gridpos( gp, x, ds );
      assert( gp[0].idx == 0 );
      interpweights( itw, gp[0] );
      interp( ext_mat_mono, itw, ext_matArray, gp[0] );
      ext_mat  = ext_mat_mono;
      ext_mat += ext_matArray[gp[0].idx];
      //
      {
        Index extmat_case = 0;
        ext2trans( incT, extmat_case, ext_mat, ds/2 );
      }
      //
      mult( evol_op, evol_opArray[gp[0].idx], incT );
      interp( abs_vec_mono, itw, abs_vecArray,gp[0] );
      temperature = interp( itw, tArray, gp[0] );
      interp( pnd_vec, itw, pnd_vecArray, gp[0] );
      for( Index i=0; i<2; i++ )
        {
          rte_pos[i] = interp( itw, ppath_step.pos(Range(ip-1,2),i), gp[0] );
          rte_los[i] = interp( itw, ppath_step.los(Range(ip-1,2),i), gp[0] );
        }
      rte_pos[2] = interp( itw, ppath_step.pos(Range(ip-1,2),2), gp[0] );
    }

  assert(isfinite(g));

  // A dirty trick to avoid copying ppath_step
  const Index np = ip+1;
  ppath_step.np  = np;

}



//! opt_propCalc
/*!
Returns the total monochromatic extinction matrix and absorption vector over all
scattering elements at a specific atmospheric location.

   \return ext_mat_mono   Output: total monochromatic extinction matrix
   \return abs_vec_mono   Output: total monochromatic absorption vector
   \param za              zenith angle of propagation direction
   \param aa              azimuthal angle of propagation
   \param scat_data_mono  as the WSV
   \param stokes_dim      as the WSV
   \param pnd_vec         vector of particle number densities (one element per
                          scattering element)
   \param rtp_temperature as the WSV

   \author Cory Davis
   \date   2004-7-16
*/
void opt_propCalc(
                  MatrixView      ext_mat_mono,
                  VectorView      abs_vec_mono,
                  const Numeric   za,
                  const Numeric   aa,
                  const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                  const Index     stokes_dim,
                  ConstVectorView pnd_vec,
                  const Numeric   rtp_temperature,
                  const Verbosity& verbosity
                  )
{
  assert( stokes_dim>=1  &&  stokes_dim<=4 );
  assert( ext_mat_mono.nrows() == stokes_dim );
  assert( ext_mat_mono.ncols() == stokes_dim );
  assert( abs_vec_mono.nelem() == stokes_dim );

  const Index N_hm = scat_data_mono.nelem();

  Matrix ext_mat_mono_spt(stokes_dim,stokes_dim);
  Vector abs_vec_mono_spt(stokes_dim);

  ext_mat_mono=0.0;
  abs_vec_mono=0.0;  

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_hm; i_ss++)
  {
      const Index N_se = scat_data_mono[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          if (pnd_vec[i_se]>0)
          {
              Index nT=scat_data_mono[i_ss][i_se].T_grid.nelem();
              if( nT > 1 )
              {
                  //set extrapolfax explicitly. should correspond to the one assumed
                  //by gridpos call in opt_propExtract.
                  Numeric extrapolfac=0.5;
                  Numeric lowlim = scat_data_mono[i_ss][i_se].T_grid[0]-
                  extrapolfac*(scat_data_mono[i_ss][i_se].T_grid[1]
                               -scat_data_mono[i_ss][i_se].T_grid[0]);
                  Numeric uplim = scat_data_mono[i_ss][i_se].T_grid[nT-1]+
                  extrapolfac*(scat_data_mono[i_ss][i_se].T_grid[nT-1]
                               -scat_data_mono[i_ss][i_se].T_grid[nT-2]);

                  if( rtp_temperature<lowlim || rtp_temperature>uplim )
                  {
                      ostringstream os;
                      os << "Atmospheric temperature (" << rtp_temperature
                      << "K) out of valid temperature\n"
                      << "range of particle optical properties ("
                      << lowlim << "-" << uplim
                      << "K) of scattering element #" << i_se << ".\n";
                      throw runtime_error( os.str() );
                  }
              }
              opt_propExtract( ext_mat_mono_spt, abs_vec_mono_spt,
                              scat_data_mono[i_ss][i_se], za, aa,
                              rtp_temperature, stokes_dim, verbosity);

              ext_mat_mono_spt *= pnd_vec[i_se];
              abs_vec_mono_spt *= pnd_vec[i_se];
              ext_mat_mono     += ext_mat_mono_spt;
              abs_vec_mono     += abs_vec_mono_spt;
          }
      }
  }
}


void opt_propExtract(
                     MatrixView     ext_mat_mono_spt,
                     VectorView     abs_vec_mono_spt,
                     const SingleScatteringData& scat_data_single,
                     const Numeric  za,
                     const Numeric  aa _U_, // avoid warning until we use ptype=10
                     const Numeric  rtp_temperature,
                     const Index    stokes_dim,
                     const Verbosity& verbosity
                     )
{
  // Temperature grid position
  GridPos t_gp;
  if( scat_data_single.T_grid.nelem() > 1)
    { gridpos( t_gp, scat_data_single.T_grid, rtp_temperature ); }


  switch (scat_data_single.ptype){

  case PTYPE_GENERAL:
    {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         as long as only PTYPE_MACROS_ISO and PTYPE_HORIZ_AL, but will be a
         problem for PTYPE_GENERAL, ie needs to be fixed BEFORE adding
         PTYPE_GENERAL support (see AUG appendix for more info).
      */

      // This is only included to remove warnings about unused variables 
      // during compilation
      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
      break;
    }
  case PTYPE_MACROS_ISO:
    {
      assert (scat_data_single.ext_mat_data.ncols() == 1);
      
      Vector itw(2);
      //interpolate over temperature
      if( scat_data_single.T_grid.nelem() > 1)
      {
        interpweights(itw, t_gp);
      }
      // In the case of randomly oriented particles the extinction matrix is 
      // diagonal. The value of each element of the diagonal is the
      // extinction cross section, which is stored in the database.
     
      ext_mat_mono_spt = 0.;
      abs_vec_mono_spt = 0;
      
      if( scat_data_single.T_grid.nelem() == 1)
        {
          ext_mat_mono_spt(0,0) = scat_data_single.ext_mat_data(0,0,0,0,0);
          abs_vec_mono_spt[0] = scat_data_single.abs_vec_data(0,0,0,0,0);
        }
      // Temperature interpolation
      else
        {
          ext_mat_mono_spt(0,0) = interp(itw,scat_data_single.ext_mat_data(0,joker,0,0,0),t_gp);
          abs_vec_mono_spt[0] = interp(itw,scat_data_single.abs_vec_data(0,joker,0,0,0),t_gp);
        }
      
      if( stokes_dim == 1 ){
        break;
      }
      
      ext_mat_mono_spt(1,1) = ext_mat_mono_spt(0,0);
      
      if( stokes_dim == 2 ){
        break;
      }
      
      ext_mat_mono_spt(2,2) = ext_mat_mono_spt(0,0);
      
      if( stokes_dim == 3 ){
        break;
      }
      
      ext_mat_mono_spt(3,3) = ext_mat_mono_spt(0,0);
      break;
    }

  case PTYPE_HORIZ_AL:
    {
      assert (scat_data_single.ext_mat_data.ncols() == 3);
      
      // In the case of horizontally oriented particles the extinction matrix
      // has only 3 independent non-zero elements ext_mat_monojj, K12=K21, and
      // K34=-K43. These values are dependent on the zenith angle of
      // propagation. The data storage format also makes use of the fact that
      // in this case K(za_sca)=K(180-za_sca).

      // 1st interpolate data by za_sca
      GridPos za_gp;
      Numeric Kjj;
      Numeric K12;
      Numeric K34;
      

      // This fix for za=90 copied from  ext_matTransform in optproperties.cc
      // (121113, PE).
      ConstVectorView this_za_datagrid = 
              scat_data_single.za_grid[ Range( 0, scat_data_single.ext_mat_data.npages() ) ];
      if( za > 90 )
        { gridpos( za_gp, this_za_datagrid, 180-za ); }
      else
        { gridpos( za_gp, this_za_datagrid, za ); }


      ext_mat_mono_spt = 0.0;
      abs_vec_mono_spt = 0.0;

      Vector itw;  

      if( scat_data_single.T_grid.nelem() == 1 )
        {
          itw.resize(2);
          interpweights( itw, za_gp );
          Kjj = interp(itw,scat_data_single.ext_mat_data(0,0,joker,0,0),za_gp);
          abs_vec_mono_spt[0]   = interp( itw, 
                          scat_data_single.abs_vec_data(0,0,joker,0,0),za_gp);
        }
      else
        {
          itw.resize(4);
          interpweights( itw, t_gp, za_gp );
          Kjj = interp(itw,scat_data_single.ext_mat_data(0,joker,joker,0,0),t_gp,za_gp);
          abs_vec_mono_spt[0]   = interp( itw, 
                       scat_data_single.abs_vec_data(0,joker,joker,0,0), t_gp,za_gp );
        }
      ext_mat_mono_spt(0,0) = Kjj;



      if( stokes_dim == 1 )
        { break; }
      
      if( scat_data_single.T_grid.nelem() == 1 )
        {
          K12=interp(itw,scat_data_single.ext_mat_data(0,0,joker,0,1),za_gp);
          abs_vec_mono_spt[1] = interp(itw,
                         scat_data_single.abs_vec_data(0,0,joker,0,1),za_gp);
        }
      else
        {
          K12=interp(itw,scat_data_single.ext_mat_data(0,joker,joker,0,1),t_gp,za_gp);
          abs_vec_mono_spt[1] = interp(itw,
                         scat_data_single.abs_vec_data(0,joker,joker,0,1),t_gp,za_gp);
        }
      ext_mat_mono_spt(1,1)=Kjj;
      ext_mat_mono_spt(0,1)=K12;
      ext_mat_mono_spt(1,0)=K12;

      if( stokes_dim == 2 ){
        break;
      }
      
      ext_mat_mono_spt(2,2)=Kjj;
      
      if( stokes_dim == 3 ){
        break;
      }
      
      if( scat_data_single.T_grid.nelem() == 1 )
        K34=interp(itw,scat_data_single.ext_mat_data(0,0,joker,0,2),za_gp);
      else
        K34=interp(itw,scat_data_single.ext_mat_data(0,joker,joker,0,2),t_gp,za_gp);
      ext_mat_mono_spt(2,3)=K34;
      ext_mat_mono_spt(3,2)=-K34;
      ext_mat_mono_spt(3,3)=Kjj;
      break;

    }
  default:
    {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
    
  }


}


//! pha_mat_singleCalc
/*!
 Returns the total phase matrix for given incident and scattered directions
. It requires a vector of particle number densities to be precalculated

 \param[out] Z               Output: phase matrix
 \param[in]  za_sca          scattered 
 \param[in]  aa_sca          and
 \param[in]  za_inc          incident
 \param[in]  aa_inc          directions
 \param[in]  scat_data_mono  as the WSV
 \param[in]  stokes_dim      as the WSV
 \param[in]  pnd_vec         vector of particle number densities at the point 
                             in question
 \param[in]  rtp_temperature as the WSV
 \author Cory Davis
 \date   2003-11-27
*/
void pha_mat_singleCalc(
                        MatrixView Z,                  
                        const Numeric    za_sca, 
                        const Numeric    aa_sca, 
                        const Numeric    za_inc, 
                        const Numeric    aa_inc,
                        const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                        const Index      stokes_dim,
                        ConstVectorView  pnd_vec,
                        const Numeric    rtp_temperature,
                        const Verbosity& verbosity
                        )
{
  assert(aa_inc>=-180 && aa_inc<=180);
  assert(aa_sca>=-180 && aa_sca<=180);

  Matrix Z_spt(stokes_dim, stokes_dim, 0.);
  Z=0.0;
  Index i_se_pnd = -1;
  // this is a loop over the different scattering elements
  for (Index i_ss = 0; i_ss<scat_data_mono.nelem(); i_ss++)
  {
      const Index N_se = scat_data_mono[i_ss].nelem();

      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          i_se_pnd++;
          if (pnd_vec[i_se_pnd]>0)
          {
              pha_mat_singleExtract(Z_spt,scat_data_mono[i_ss][i_se],za_sca,aa_sca,za_inc,
                                    aa_inc,rtp_temperature,stokes_dim,verbosity);
              Z_spt*=pnd_vec[i_se];
              Z+=Z_spt;
          }
      }
  }
}



//! Extract the phase matrix from a monochromatic SingleScatteringData object
/*!
  Given a monochromatic SingleScatteringData object, incident and 
  scattered directions, and the temperature, this function returns the phase
  matrix in the laboratory frame
  \param[out] Z_spt the phase matrix
  \param[in]  scat_data_single a monochromatic SingleScatteringData object
  \param[in]  za_sca 
  \param[in]  aa_sca
  \param[in]  za_inc
  \param[in]  aa_inc
  \param[in]  rtp_temperature
  \param[in]  stokes_dim

  \author Cory Davis
  \date 2004-07-16

*/
void pha_mat_singleExtract(
                           MatrixView Z_spt,
                           const SingleScatteringData& scat_data_single,
                           const Numeric za_sca,
                           const Numeric aa_sca,
                           const Numeric za_inc,
                           const Numeric aa_inc,
                           const Numeric rtp_temperature,
                           const Index   stokes_dim,
                           const Verbosity& verbosity
                           )                       
{
  switch (scat_data_single.ptype){

    case PTYPE_GENERAL:
    {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         as long as only PTYPE_MACROS_ISO and PTYPE_HORIZ_AL, but will be a
         problem for PTYPE_GENERAL, ie needs to be fixed BEFORE adding
         PTYPE_GENERAL support (see AUG appendix for more info).
      */

      // to remove warnings during compilation.
      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
      break;
    }
  case PTYPE_MACROS_ISO:
    {
      // Calculate the scattering and interpolate the data on the scattering
      // angle:
      
      Vector pha_mat_int(6);
      Numeric theta_rad;
            
      // Interpolation of the data on the scattering angle:
      interp_scat_angle_temperature(pha_mat_int, theta_rad, 
                                    scat_data_single, za_sca, aa_sca, 
                                    za_inc, aa_inc, rtp_temperature);
      
      // Caclulate the phase matrix in the laboratory frame:
      pha_mat_labCalc(Z_spt, pha_mat_int, za_sca, aa_sca, za_inc, aa_inc, theta_rad);
      
      break;
    }

  case PTYPE_HORIZ_AL:
    //Data is already stored in the laboratory frame, but it is compressed
    //a little.  Details elsewhere
    {
      assert (scat_data_single.pha_mat_data.ncols()==16);
      Numeric delta_aa=aa_sca-aa_inc+(aa_sca-aa_inc<-180)*360-
        (aa_sca-aa_inc>180)*360;//delta_aa corresponds to the "pages" 
                                //dimension of pha_mat_data
      GridPos t_gp;
      GridPos za_sca_gp;
      GridPos delta_aa_gp;
      GridPos za_inc_gp;
      Vector itw;
      
      gridpos(delta_aa_gp,scat_data_single.aa_grid,abs(delta_aa));
      if (za_inc>90)
        {
          gridpos(za_inc_gp,scat_data_single.za_grid,180-za_inc);
          gridpos(za_sca_gp,scat_data_single.za_grid,180-za_sca);
        }
      else
        {
          gridpos(za_inc_gp,scat_data_single.za_grid,za_inc);
          gridpos(za_sca_gp,scat_data_single.za_grid,za_sca);
        }

      if( scat_data_single.T_grid.nelem() == 1 )
        {
          itw.resize(8);
          interpweights(itw,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(0,0)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,0),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      else
        {
          itw.resize(16);
          gridpos(t_gp,scat_data_single.T_grid,rtp_temperature);
          interpweights(itw,t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(0,0)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,0),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      
      if( stokes_dim == 1 ){
        break;
      }

      if( scat_data_single.T_grid.nelem() == 1 )
        {
          Z_spt(0,1)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,1),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,0)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,4),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,1)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,5),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      else
        {
          Z_spt(0,1)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,1),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,0)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,4),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,1)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,5),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }

      if( stokes_dim == 2 ){
        break;
      }

      if ((za_inc<=90 && delta_aa>=0)||(za_inc>90 && delta_aa<0))
        {
          if( scat_data_single.T_grid.nelem() == 1 )
            {
              Z_spt(0,2)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,2),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,2)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,6),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,0)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,8),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,1)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,9),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
            }
          else
            {
              Z_spt(0,2)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,2),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,2)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,6),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,0)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,8),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,1)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,9),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
            }
        }
      else
        {
          if( scat_data_single.T_grid.nelem() == 1 )
            {
              Z_spt(0,2)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,2),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,2)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,6),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,0)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,8),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,1)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,9),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
            }
          else
            {
              Z_spt(0,2)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,2),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,2)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,6),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,0)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,8),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(2,1)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,9),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
            }
        }                             
      if( scat_data_single.T_grid.nelem() == 1 )
        Z_spt(2,2)=interp(itw,scat_data_single.pha_mat_data
                          (0,0,joker,joker,joker,0,10),
                          za_sca_gp,delta_aa_gp,za_inc_gp);
      else
        Z_spt(2,2)=interp(itw,scat_data_single.pha_mat_data
                          (0,joker,joker,joker,joker,0,10),
                          t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);

      if( stokes_dim == 3 ){
        break;
      }

      if ((za_inc<=90 && delta_aa>=0)||(za_inc>90 && delta_aa<0))
        {
          if( scat_data_single.T_grid.nelem() == 1 )
            {
              Z_spt(0,3)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,3),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,3)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,7),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,0)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,12),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,1)=interp(itw,scat_data_single.pha_mat_data
                                (0,0,joker,joker,joker,0,13),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
            }
          else
            {
              Z_spt(0,3)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,3),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,3)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,7),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,0)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,12),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,1)=interp(itw,scat_data_single.pha_mat_data
                                (0,joker,joker,joker,joker,0,13),
                                t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
            }
        }
      else
        {
          if( scat_data_single.T_grid.nelem() == 1 )
            {
              Z_spt(0,3)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,3),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,3)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,7),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,0)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,12),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,1)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,0,joker,joker,joker,0,13),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
            }
          else
            {
              Z_spt(0,3)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,3),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(1,3)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,7),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,0)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,12),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
              Z_spt(3,1)=-interp(itw,scat_data_single.pha_mat_data
                                 (0,joker,joker,joker,joker,0,13),
                                 t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
            }
        }

      if( scat_data_single.T_grid.nelem() == 1 )
        {
          Z_spt(2,3)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,11),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,2)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,14),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,3)=interp(itw,scat_data_single.pha_mat_data
                            (0,0,joker,joker,joker,0,15),
                            za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      else
        {
          Z_spt(2,3)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,11),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,2)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,14),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,3)=interp(itw,scat_data_single.pha_mat_data
                            (0,joker,joker,joker,joker,0,15),
                            t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      break;
      
    }  
  default:
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    
  }
}




//! Sample_los
/*!
  FIXME: 2011-06-17 Documentation removed by Gerrit (severely out of date)

   \param[out]    new_rte_los     incident line of sight for subsequent 
   \param[out]    g_los_csc_theta probability density for the chosen
                                  direction multiplied by sin(za)
   \param[out]    Z
   \param[in,out] rng             Rng random number generator instance
   \param[in]     rte_los         incident line of sight for subsequent 
                                  ray-tracing.                     
   \param[in]     scat_data_mono
   \param[in]     stokes_dim
   \param[in]     pnd_vec
   \param[in]     anyptype30
   \param[in]     Z11maxvector
   \param[in]     Csca
   \param[in]     rtp_temperature

   \author Cory Davis
   \date   2003-06-19
*/

void Sample_los (
                 VectorView       new_rte_los,
                 Numeric&         g_los_csc_theta,
                 MatrixView       Z,
                 Rng&             rng,
                 ConstVectorView  rte_los,
                 const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                 const Index      stokes_dim,
                 ConstVectorView  pnd_vec,
                 const bool       anyptype30,
                 ConstVectorView  Z11maxvector,
                 const Numeric    Csca,
                 const Numeric    rtp_temperature,
                 const Verbosity& verbosity
                 )
{
  Numeric Z11max=0;
  bool tryagain=true;

  Vector sca_dir;
  mirror_los( sca_dir, rte_los, 3 );
      
  // Rejection method http://en.wikipedia.org/wiki/Rejection_sampling
  if(anyptype30)
    {
      Index np=pnd_vec.nelem();
      assert(TotalNumberOfElements(scat_data_mono)==np);
      for(Index i=0;i<np;i++)
        {
          Z11max+=Z11maxvector[i]*pnd_vec[i];
        }
    }
  else
    {
      Matrix dummyZ(stokes_dim,stokes_dim);
      //The following is based on the assumption that the maximum value of the 
      //phase matrix for a given scattered direction is for forward scattering
      pha_mat_singleCalc( dummyZ, sca_dir[0], sca_dir[1], sca_dir[0], sca_dir[1],
                          scat_data_mono, stokes_dim, pnd_vec, rtp_temperature,
                          verbosity);
      Z11max=dummyZ(0,0);
    }  
  ///////////////////////////////////////////////////////////////////////  
  while(tryagain)
    {
      new_rte_los[0] = acos(1-2*rng.draw())*RAD2DEG;
      new_rte_los[1] = rng.draw()*360-180;

      //Calculate Phase matrix////////////////////////////////
      Vector inc_dir;
      mirror_los( inc_dir, new_rte_los, 3 );
      
      pha_mat_singleCalc( Z, sca_dir[0], sca_dir[1], inc_dir[0], inc_dir[1],
                          scat_data_mono, stokes_dim, pnd_vec, rtp_temperature,
                          verbosity );
      
      if (rng.draw()<=Z(0,0)/Z11max)//then new los is accepted
        {
          tryagain=false;
        }
    }
  g_los_csc_theta =Z(0,0)/Csca;
}


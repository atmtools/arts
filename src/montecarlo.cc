/* copyright (C) 2003-2007 Cory Davis <cory@met.ed.ac.uk>
                            
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

  \brief  functions used by ScatteringMonteCarlo
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "auto_md.h"
#include "montecarlo.h"

#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include "sstream.h"
#endif

//! clear_rt_vars_at_gp
/*! 

  Calculates a bunch of atmospheric variables at the end of a ppath.
   
\author Cory Davis
\date 2005-02-19?
*/  


void clear_rt_vars_at_gp(
                         MatrixView&             ext_mat_mono,
                         VectorView&             abs_vec_mono,
                         Numeric&                temperature,
                         const Agenda&           opt_prop_gas_agenda,
                         const Agenda&           abs_scalar_gas_agenda,
                         const Index&            f_index,
                         const GridPos&          gp_p,
                         const GridPos&          gp_lat,
                         const GridPos&          gp_lon,
                         const ConstVectorView   p_grid,
                         const ConstVectorView   lat_grid,
                         const ConstVectorView   lon_grid,
                         const ConstTensor3View  t_field,
                         const ConstTensor4View  vmr_field)
{
  const Index   ns = vmr_field.nbooks();
  Vector t_vec(1);  //vectors are required by interp_atmfield_gp2itw etc.
  Vector   p_vec(1);//may not be efficient with unecessary vectors
  Matrix itw_p(1,2);
  ArrayOfGridPos ao_gp_p(1),ao_gp_lat(1),ao_gp_lon(1);
  Matrix   vmr_mat(ns,1), itw_field;
  
  //local versions of workspace variables
  Matrix local_abs_scalar_gas,local_abs_vec;
  Tensor3 local_ext_mat;
  ao_gp_p[0]=gp_p;
  ao_gp_lat[0]=gp_lat;
  ao_gp_lon[0]=gp_lon;
  
  // Determine the pressure 

  interpweights( itw_p, ao_gp_p );
  itw2p( p_vec, p_grid, ao_gp_p, itw_p );
  
  // Determine the atmospheric temperature and species VMR 


  interp_atmfield_gp2itw( itw_field, 3, p_grid, 
                          lat_grid, lon_grid, 
                          ao_gp_p, ao_gp_lat, ao_gp_lon );
  //
  interp_atmfield_by_itw( t_vec,  3, p_grid, 
                          lat_grid, lon_grid, t_field, 
                          "t_field", ao_gp_p, ao_gp_lat, ao_gp_lon, 
                          itw_field );
  // 
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( vmr_mat(is, joker), 3, p_grid, 
                              lat_grid, lon_grid, 
                              vmr_field(is, joker, joker, joker), 
                              "vmr_field", ao_gp_p, ao_gp_lat, 
                              ao_gp_lon, itw_field );
    }


  temperature = t_vec[0];
  
  //calcualte absorption coefficient
  abs_scalar_gas_agendaExecute( local_abs_scalar_gas,f_index,p_vec[0],
    temperature,vmr_mat(joker,0),abs_scalar_gas_agenda,true );
  opt_prop_gas_agendaExecute( local_ext_mat, local_abs_vec, local_abs_scalar_gas,
                              opt_prop_gas_agenda, true );
  ext_mat_mono=local_ext_mat(0, Range(joker), Range(joker));
  abs_vec_mono=local_abs_vec(0,Range(joker));
}



//! cloudy_rt_vars_at_gp
/*! 

  Calculates a bunch of atmospheric variables at the end of a ppath.
   
\author Cory Davis
\date 2005-02-19?
*/  


void cloudy_rt_vars_at_gp(
                          MatrixView&           ext_mat_mono,
                          VectorView&           abs_vec_mono,
                          VectorView&           pnd_vec,
                          Numeric&              temperature,
                          const Agenda&         opt_prop_gas_agenda,
                          const Agenda&         abs_scalar_gas_agenda,
                          const Index&          stokes_dim,
                          const Index&          f_index,
                          const GridPos         gp_p,
                          const GridPos         gp_lat,
                          const GridPos         gp_lon,
                          const ConstVectorView    p_grid_cloud,
                          const ConstVectorView    lat_grid_cloud,
                          const ConstVectorView    lon_grid_cloud,
                          const ConstTensor3View   t_field_cloud,
                          const ConstTensor4View   vmr_field_cloud,
                          const Tensor4&           pnd_field,
                          const ArrayOfSingleScatteringData& scat_data_mono,
                          const ArrayOfIndex&                cloudbox_limits,
                          const Vector&             rte_los
                          )

{
  const Index   ns = vmr_field_cloud.nbooks();
  const Index N_pt = pnd_field.nbooks();
  Matrix  pnd_ppath(N_pt,1);
  Vector t_ppath(1);
  Vector   p_ppath(1);//may not be efficient with unecessary vectors
  Matrix itw_p(1,2);
  ArrayOfGridPos ao_gp_p(1),ao_gp_lat(1),ao_gp_lon(1);
  Matrix   vmr_ppath(ns,1), itw_field;
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);
  Vector abs_vec_part(stokes_dim, 0.0);
  Numeric scat_za,scat_aa;

  //local versions of workspace variables
  Matrix local_abs_scalar_gas,local_abs_vec;
  Tensor3 local_ext_mat;
  //Numeric local_rte_pressure;
  //Vector local_rte_vmrlist;

  ao_gp_p[0]=gp_p;
  ao_gp_lat[0]=gp_lat;
  ao_gp_lon[0]=gp_lon;
  

  cloud_atm_vars_by_gp(p_ppath,t_ppath,vmr_ppath,pnd_ppath,ao_gp_p,
                       ao_gp_lat,ao_gp_lon,cloudbox_limits,p_grid_cloud,
                       lat_grid_cloud,lon_grid_cloud,t_field_cloud,
                       vmr_field_cloud,pnd_field);
   
  //local_rte_pressure    = p_ppath[0];
  temperature = t_ppath[0];
  //rte_vmr_list    = vmr_ppath(joker,0);
  abs_scalar_gas_agendaExecute( local_abs_scalar_gas,f_index,p_ppath[0],
        temperature,vmr_ppath(joker,0),abs_scalar_gas_agenda,true );
  opt_prop_gas_agendaExecute( local_ext_mat, local_abs_vec, local_abs_scalar_gas,
                              opt_prop_gas_agenda, true );
  ext_mat_mono=local_ext_mat(0, Range(joker), Range(joker));
  abs_vec_mono=local_abs_vec(0,Range(joker));
  ext_mat_part=0.0;
  abs_vec_part=0.0;
  scat_za=180-rte_los[0];
  scat_aa=rte_los[1]+180;
  //Make sure scat_aa is between -180 and 180
  if (scat_aa>180){scat_aa-=360;}
  //
  //opt_prop_part_agenda.execute( true );
  //use pnd_ppath and ext_mat_spt to get extmat (and similar for abs_vec
  pnd_vec=pnd_ppath(joker, 0);
  opt_propCalc(ext_mat_part,abs_vec_part,scat_za,scat_aa,scat_data_mono,
               stokes_dim, pnd_vec, temperature);
  
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
  \param vmr          Output: a n_species by n matrix of VMRs
  \param pnd          Output: a n_ptypes by n matrix of VMRs
  \param gp_p         an array of pressre gridpoints
  \param gp_lat       an array of latitude gridpoints
  \param gp_lon       an array of longitude gridpoints
  \param cloudbox_limits  the WSV
  \param p_grid_cloud the subset of the p_grid corresponding to the cloudbox
  \param lat_grid_cloud the subset of the lat_grid corresponding to the cloudbox
  \param lon_grid_cloud the subset of the lon_grid corresponding to the cloudbox
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
                          const ArrayOfIndex& cloudbox_limits,
                          const ConstVectorView p_grid_cloud,
                          const ConstVectorView lat_grid_cloud,
                          const ConstVectorView lon_grid_cloud,
                          const ConstTensor3View   t_field_cloud,
                          const ConstTensor4View   vmr_field_cloud,
                          const ConstTensor4View   pnd_field
                          )
{
  Index np=gp_p.nelem();
  assert(pressure.nelem()==np);
  Index ns=vmr_field_cloud.nbooks();
  Index N_pt=pnd_field.nbooks();
  ArrayOfGridPos gp_p_cloud = gp_p;
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
  
  fix_gridpos_at_boundary(gp_p_cloud, p_grid_cloud.nelem()); 
  fix_gridpos_at_boundary(gp_lat_cloud, lat_grid_cloud.nelem()); 
  fix_gridpos_at_boundary(gp_lon_cloud, lon_grid_cloud.nelem());
  
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
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid_cloud, 
                          lat_grid_cloud, lon_grid_cloud, 
                          gp_p_cloud, gp_lat_cloud, gp_lon_cloud );
  //
  interp_atmfield_by_itw( temperature,  atmosphere_dim, p_grid_cloud, 
                          lat_grid_cloud, lon_grid_cloud, t_field_cloud, 
                          "t_field", gp_p_cloud, gp_lat_cloud, gp_lon_cloud, 
                          itw_field );
  // 
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( vmr(is, joker), atmosphere_dim, p_grid_cloud, 
                              lat_grid_cloud, lon_grid_cloud, 
                              vmr_field_cloud(is, joker, joker, joker), 
                              "vmr_field", gp_p_cloud, gp_lat_cloud, 
                              gp_lon_cloud, itw_field );
    }
  
  //Determine the particle number density for every particle type at 
  // each propagation path point
  for( Index ip=0; ip<N_pt; ip++ )
    {
      // if grid positions still outside the range the propagation path step 
      // must be outside the cloudbox and pnd is set to zero
      interp_atmfield_by_itw( pnd(ip, joker), atmosphere_dim,
                              p_grid_cloud, lat_grid_cloud, lon_grid_cloud, 
                              pnd_field(ip, joker, joker,  joker), 
                              "pnd_field", gp_p_cloud, gp_lat_cloud, 
                              gp_lon_cloud, itw_field );
    }
}


//! Cloudbox_ppath_calc
/*! 
  This function performs the same task as ppath_calc, except inside the
  cloudbox.  It has been derived from the clear sky version.  See the 
  online help (arts -d FUNCTION_NAME) for a description of parameters.
   
\author Cory Davis
\date 2003-06-19
*/  

void Cloudbox_ppathCalc(
                        //  Output:
                        Ppath&          ppath,
                        Ppath&          ppath_step,
                        //  Input:
                        const Agenda&         ppath_step_agenda,
                        const Index&          atmosphere_dim,
                        const Vector&         p_grid,
                        const Vector&         lat_grid,
                        const Vector&         lon_grid,
                        const Tensor3&        z_field,
                        const Matrix&         r_geoid,
                        const Matrix&         z_surface,
                        const ArrayOfIndex&   cloudbox_limits,
                        const Vector&         rte_pos,
                        const Vector&         rte_los,
                        const Index& z_field_is_1D
)
{
  // This function is a WSM but it is normally only called from RteCalc. 
  // For that reason, this function does not repeat input checks that are
  // performed in RteCalc, it only performs checks regarding the sensor 
  // position and LOS.

  //--- Check input -----------------------------------------------------------

  // Sensor position and LOS
  //
  chk_vector_length( "rte_pos", rte_pos, atmosphere_dim );
  chk_if_over_0( "sensor radius", rte_pos[0] );
  if( atmosphere_dim < 3 )
    {
        ostringstream os;
        os << "cloudbox_ppath_calc only works for a 3D atmosphere";
        throw runtime_error( os.str() );
    }
  else
    {
      chk_if_in_range( "sensor latitude", rte_pos[1], -90., 90. );
      chk_if_in_range( "sensor longitude", rte_pos[2], -360., 360. );
      chk_vector_length( "rte_los", rte_los, 2 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], 0., 180. );
      chk_if_in_range( "sensor azimuth angle", rte_los[1], -180., 180. );
    }
  
  //--- End: Check input ------------------------------------------------------


  // Some messages
  out2 << "  -------------------------------------\n";
  out2 << "  sensor radius          : " << rte_pos[0]/1e3 << " km\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor longitude       : " << rte_pos[2] << "\n";
  out2 << "  sensor zenith angle    : " << rte_los[0] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor azimuth angle   : " << rte_los[1] << "\n";
  
  
  
  // Initiate the partial Ppath structure. 

  //
  cloudbox_ppath_start_stepping( ppath_step, atmosphere_dim, p_grid, lat_grid, 
                                 lon_grid, z_field, r_geoid, z_surface, rte_pos, rte_los, 
                                 z_field_is_1D );

  out2 << "  -------------------------------------\n";

  // Perform propagation path steps until the starting point is found, which
  // is flagged by ppath_step by setting the background field.
  //
  // The results of each step, returned by ppath_step_agenda as a new 
  // ppath_step, are stored as an array of Ppath structures.
  //
  Array<Ppath>   ppath_array;
  ppath_array.push_back( ppath_step );
  // 
  Index   np = ppath_step.np;   // Counter for number of points of the path
  Index   istep = 0;            // Counter for number of steps
//

  
  while( !ppath_what_background( ppath_step ) )
    {

      // Call ppath_step agenda. 
      // The new path step is added to *ppath_array* last in the while block
      //
      istep++;
      //
      ppath_step_agendaExecute(ppath_step, atmosphere_dim, p_grid,
                               lat_grid, lon_grid, z_field, r_geoid, z_surface,
                               ppath_step_agenda, true);

      // Before everything is tested carefully, we consider more than 1000
      // path points to be an indication on that the calcululations have
      // got stuck in an infinite loop.
      if( istep > 5000 )
        {
          throw runtime_error(
             "5000 path points have been reached. Is this an infinite loop?" );
        }
      //     cout << "istep = " << istep << "\n";
      // Number of points in returned path step
      const Index n = ppath_step.np;

      // Increase the total number
      np += n - 1;
      
      // Check if there is an intersection with the cloud box boundary,
        //remembering that this function will only be called within the
        //cloud box
      
      double ipos = double( ppath_step.gp_p[n-1].idx ) + 
                                                    ppath_step.gp_p[n-1].fd[0];
      if( ipos <= double( cloudbox_limits[0] )  || 
                                        ipos >= double( cloudbox_limits[1] ) )
        { ppath_set_background( ppath_step, 3 ); }
      else   
          {
          ipos = double( ppath_step.gp_lat[n-1].idx ) + 
                                              ppath_step.gp_lat[n-1].fd[0];
          if( ipos <= double( cloudbox_limits[2] )  || 
                                    ipos >= double( cloudbox_limits[3] ) )
             { ppath_set_background( ppath_step, 3 ); }
          else
         {
           ipos = double( ppath_step.gp_lon[n-1].idx ) + 
                                              ppath_step.gp_lon[n-1].fd[0];
              if( ipos <= double( cloudbox_limits[4] )  || 
                                    ipos >= double( cloudbox_limits[5] ) )
                    { ppath_set_background( ppath_step, 3 ); } 
            }
        }
        
      // Put new ppath_step in ppath_array
      ppath_array.push_back( ppath_step );
      //     cout <<"what background? "<< ppath_what_background( ppath_step )<< "\n";
    } // End path steps
  
 
  // Combine all structures in ppath_array to form the return Ppath structure.
  //
  ppath_init_structure( ppath, atmosphere_dim, np );
  //
  np = 0;   // Now used as counter for points moved to ppath
  //
  for( Index i=0; i<ppath_array.nelem(); i++ )
    {
      // For the first structure, the first point shall be included, but the
      // first structure can also be empty. 
      // For later structures, the first point shall not be included, but
      // there will always be at least two points.
      // Only the first structure can be empty.

      Index n = ppath_array[i].np;

      if( n )
        {
          // First index to include
          Index i1 = 1;
          if( i == 0 )
            { i1 = 0; }
          else
            { assert( n > 1 ); }

          // Vectors and matrices that can be handled by ranges.
          ppath.z[ Range(np,n-i1) ] = ppath_array[i].z[ Range(i1,n-i1) ];
          ppath.pos( Range(np,n-i1), joker ) = 
                                   ppath_array[i].pos( Range(i1,n-i1), joker );
          ppath.los( Range(np,n-i1), joker ) = 
                                   ppath_array[i].los( Range(i1,n-i1), joker );

          // For i==1, there is no defined l_step. For higher i, all 
          // values in l_step shall be copied.
          if( i > 0 )
            { ppath.l_step[ Range(np-1,n-1) ] = ppath_array[i].l_step; }

          // Grid positions must be handled by a loop
          for( Index j=i1; j<n; j++ )
            { ppath.gp_p[np+j-i1] = ppath_array[i].gp_p[j]; }
          if( atmosphere_dim >= 2 )
            {
              for( Index j=i1; j<n; j++ )
                { ppath.gp_lat[np+j-i1] = ppath_array[i].gp_lat[j]; }
            }
          if( atmosphere_dim == 3 )
            {
              for( Index j=i1; j<n; j++ )
                { ppath.gp_lon[np+j-i1] = ppath_array[i].gp_lon[j]; }
            }

          // Fields just set once
          if( ppath_array[i].tan_pos.nelem() )
            {
              ppath.tan_pos.resize( ppath_array[i].tan_pos.nelem() );
              ppath.tan_pos               = ppath_array[i].tan_pos; 
            }
          if( ppath_array[i].geom_tan_pos.nelem() )
            {
              ppath.geom_tan_pos.resize( ppath_array[i].tan_pos.nelem() );
              ppath.geom_tan_pos          = ppath_array[i].geom_tan_pos; 
            }

          // Increase number of points done
          np += n - i1;
         
        }
    }  
  ppath.method     = ppath_step.method;
  ppath.refraction = ppath_step.refraction;
  ppath.constant   = ppath_step.constant;
  ppath.background = ppath_step.background;

  

  out3 << "  number of path steps  : " << istep           << "\n";
  out3 << "  number of path points : " << ppath.z.nelem() << "\n";


  // If refraction has been considered, make a simple check that the
  // refraction at the top of the atmosphere is sufficiently close to 1.
  if( ppath.refraction  &&  min( z_field(z_field.npages()-1,0,0) ) < 60e3 )
    {
      out2 << "  *** WARNING****\n" 
           << "  The calculated propagation path can be inexact as the "
           << "atmosphere\n  only extends to " 
           <<  min( z_field(z_field.npages()-1,0,0) ) << " km. \n" 
           << "  The importance of this depends on the observation "
           << "geometry.\n  It is recommended that the top of the atmosphere "
           << "is not below 60 km.\n";
    }
}






//! Cloudbox_ppath_rteCalc

/*!
   This function was written to eliminate some repeated code in 
   ScatteringMonteCarlo. Basically this function performs internal and external 
   ppath calculations, calculates incoming radiances, and provides scattering 
   properties and other derived data that is required in the monte carlo 
   simulation for a given propagation direction.


   Parameters not already described elsewhere:
   
   \param[out] ppathcloud
   \param[out] ppath
   \param[out] ppath_step
   \param[out] rte_pos
   \param[out] rte_los
   \param[out] cum_l_step
   \param[out] TArray
   \param[out] ext_matArray
   \param[out] abs_vecArray
   \param[out] t_ppath
   \param[out] ext_mat
   \param[out] abs_vec
   \param[out] rte_pressure
   \param[out] rte_temperature
   \param[out] rte_vmr_list
   \param[out] iy
   \param[out] pnd_ppath
   \param[in]  ppath_step_agenda
   \param[in]  atmosphere_dim
   \param[in]  p_grid
   \param[in]  lat_grid
   \param[in]  lon_grid
   \param[in]  z_field
   \param[in]  r_geoid
   \param[in]  z_surface
   \param[in]  cloudbox_limits
   \param[in]  opt_prop_gas_agenda
   \param[in]  abs_scalar_gas_agenda
   \param[in]  f_index
   \param[in]  stokes_dim
   \param[in]  t_field
   \param[in]  vmr_field
   \param[in]  rte_agenda
   \param[in]  iy_space_agenda
   \param[in]  surface_prop_agenda
   \param[in]  iy_cloudbox_agenda
   \param[in]  f_grid
   \param[in]  pnd_field
   \param[in]  scat_data_mono
   \param[in]  z_field_is_1D
   \param[in]  record_ppathcloud   A flag that determines whether ppaths within the 
                                 cloud box are recorded for later plotting in 
                                 MATLAB.
   \param[in]  record_ppath        Similar, but for ppaths determining incoming 
                                 radiances at the cloudbox boundaries.  For 
                                 normal use both of these parameterss should be 
                                 set to 0.
   \param[in]  photon_number       used to label ppath files
   \param[in]  scattering_order    used to label ppath files

   \author Cory Davis
   \date   2003-06-20
*/

void Cloudbox_ppath_rteCalc(
                             Ppath&                ppathcloud,
                             Ppath&                ppath,
                             Ppath&                ppath_step,
                             Vector&               rte_pos,
                             Vector&               rte_los,
                             Vector&               cum_l_step,
                             ArrayOfMatrix&        TArray,
                             ArrayOfMatrix&        ext_matArray,
                             ArrayOfVector&        abs_vecArray,
                             Vector&               t_ppath,
                             //Vector&               scat_za_grid,
                             //Vector&               scat_aa_grid,
                             Tensor3&              ext_mat,
                             Matrix&               abs_vec,
                             Numeric&              rte_pressure,
                             Numeric&              rte_temperature,
                             Vector&               rte_vmr_list,
                             Matrix&               iy,
                             //Matrix&               i_space,
                             //Matrix&               ground_emission,
                             //Matrix&               ground_los, 
                             //Tensor4&              ground_refl_coeffs,
                             Matrix&               pnd_ppath,
                             const Agenda&         ppath_step_agenda,
                             const Index&          atmosphere_dim,
                             const Vector&         p_grid,
                             const Vector&         lat_grid,
                             const Vector&         lon_grid,
                             const Tensor3&        z_field,
                             const Matrix&         r_geoid,
                             const Matrix&         z_surface,
                             const ArrayOfIndex&   cloudbox_limits,
                             const Index&          record_ppathcloud,
                             const Index&          record_ppath,
                             const Agenda&         opt_prop_gas_agenda,
                             const Agenda&         abs_scalar_gas_agenda,
                             const Index&          f_index,
                             const Index&          stokes_dim,
                             const Tensor3&        t_field,
                             const Tensor4&        vmr_field,
                             const Agenda&         rte_agenda,
                             const Agenda&         iy_space_agenda,
                             const Agenda&         surface_prop_agenda,
                             const Agenda&         iy_cloudbox_agenda,
                             const Vector&         f_grid,
                             const Index&          photon_number,
                             const Index&          scattering_order,
                             const Tensor4&        pnd_field,
                             const ArrayOfSingleScatteringData& scat_data_mono,
                             const Index& z_field_is_1D
)

{
  // Assign dummies for variables associated with sensor.
  Vector   mblock_za_grid_dummy(1);
           mblock_za_grid_dummy[0] = 0;
  Vector   mblock_aa_grid_dummy(0);
  //Index    antenna_dim_dummy = 1;
  Sparse   sensor_response_dummy;

  // Dummy for measurement vector
  Vector   y_dummy(0);

  const Index cloudbox_on_dummy=0;
  Matrix sensor_pos(1,3);
  Matrix sensor_los(1,2);
  Tensor7 scat_i_p_dummy;
  Tensor7 scat_i_lat_dummy;
  Tensor7 scat_i_lon_dummy;
 
  //  cout << "Cloudbox_ppathCalc\n";
  Cloudbox_ppathCalc(ppathcloud,ppath_step,ppath_step_agenda,atmosphere_dim,
                     p_grid,lat_grid,lon_grid,z_field,r_geoid,z_surface,
                     cloudbox_limits, rte_pos,rte_los,z_field_is_1D);
  
  //ppath_calc(ppathcloud,ppath_step,ppath_step_agenda,atmosphere_dim,
  //           p_grid,lat_grid,lon_grid,z_field,r_geoid,z_surface,1,
  //           cloudbox_limits, rte_pos,rte_los,0);

  if (record_ppathcloud)
    {
      //Record ppathcloud.  This is useful for debugging and educational 
      //purposes.  It would be completely daft to leave this on in a 
      //real calculation
      ostringstream filename;
      filename <<"ppathcloud" << photon_number <<"_"<<scattering_order;
      String longfilename;
      filename_xml(longfilename,filename.str());
      xml_write_to_file(longfilename, ppathcloud, FILE_TYPE_ASCII);
    }

  cum_l_stepCalc(cum_l_step,ppathcloud);
  Range p_range(cloudbox_limits[0], 
                cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range(cloudbox_limits[2], 
                cloudbox_limits[3]-cloudbox_limits[2]+1);
  Range lon_range(cloudbox_limits[4], 
                cloudbox_limits[5]-cloudbox_limits[4]+1);
  //Calculate array of transmittance matrices
  TArrayCalc(TArray, ext_matArray, abs_vecArray, t_ppath, ext_mat, abs_vec, 
             rte_pressure, rte_temperature, 
             rte_vmr_list, pnd_ppath, ppathcloud, opt_prop_gas_agenda, 
             abs_scalar_gas_agenda, f_index, stokes_dim, 
             p_grid[p_range], lat_grid[lat_range], lon_grid[lon_range], 
             t_field(p_range,lat_range,lon_range), 
             vmr_field(joker,p_range,lat_range,lon_range),
             pnd_field,scat_data_mono,cloudbox_limits);

  iy_calc_no_jacobian(iy, ppath, ppath_step, ppath_step_agenda, 
                      rte_agenda, iy_space_agenda, surface_prop_agenda,
                      iy_cloudbox_agenda, atmosphere_dim, p_grid, lat_grid,
                      lon_grid, z_field, t_field, vmr_field, r_geoid, z_surface,
                      cloudbox_on_dummy, cloudbox_limits,
                      rte_pos, rte_los, f_grid,stokes_dim, 1);

  for (Index i = 0;i<stokes_dim;i++){assert(!isnan(iy(0,i)));}
  
  if (record_ppath)
    {
      //Record ppathcloud.  This is useful for debugging and educational purposes.  It would
      //be completely daft to leave this on in a real calculation
      ostringstream filename;
      filename <<"ppath" << photon_number <<"_"<<scattering_order;
      String longfilename;
      filename_xml(longfilename,filename.str());
      xml_write_to_file(longfilename, ppath, FILE_TYPE_ASCII);
    }
  
}


/*
//! cloudbox_ppath_start_stepping


   This function was derived from Patrick's 'ppath_start_stepping' function.  
   It has been adapted for use within the cloud box and is intended for a 
   3D atmosphere only

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   r_geoid           The geoid radius.
   \param   z_ground          Ground altitude.
   \param   cloudbox_on       Flag to activate the cloud box.
   \param   cloudbox_limits   Index limits of the cloud box.
   \param   rte_pos             The position of the sensor.
   \param   rte_los             The line-of-sight of the sensor.

   \author Cory Davis (derived from ppath_start_stepping (Patrick Eriksson))
   \date   2003-06-19
*/
void cloudbox_ppath_start_stepping(
                                   Ppath&                ppath,
                                   const Index&          atmosphere_dim,
                                   ConstVectorView       p_grid,
                                   ConstVectorView       lat_grid,
                                   ConstVectorView       lon_grid,
                                   ConstTensor3View      z_field,
                                   ConstMatrixView       r_geoid,
                                   ConstMatrixView       z_ground _U_,
                                   ConstVectorView       rte_pos,
                                   ConstVectorView       rte_los,
                                   const Index&          z_field_is_1D)
{
  // This function contains no checks or asserts as it is only a sub-function
  // to ppathCalc where the input is checked carefully.

  // Allocate the ppath structure
  ppath_init_structure(  ppath, atmosphere_dim, 1 );

  // Number of pressure levels
  const Index np = p_grid.nelem();

  // The different atmospheric dimensionalities are handled seperately

  if( atmosphere_dim == 3 )
  {

      // Is the sensor inside the latitude and longitude ranges of the
      // model atmosphere, and below the top of the atmosphere? If
      // yes, is_inside = true. Store geoid and ground radii, grid
      // position and interpolation weights for later use.
      //
      double rv_geoid=-1;
      //double rv_ground=-1;  // -1 to avoid compiler warnings
      GridPos   gp_lat, gp_lon;
      Vector    itw(4);
      
      
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      gridpos( gp_lon, lon_grid, rte_pos[2] );
      interpweights( itw, gp_lat, gp_lon );
      
      rv_geoid  = interp( itw, r_geoid, gp_lat, gp_lon );
      //rv_ground = rv_geoid + interp( itw, z_ground, gp_lat, gp_lon );

//          out2 << "  sensor altitude        : " << (rte_pos[0]-rv_geoid)/1e3 
//               << " km\n";
//
      // If downwards, calculate geometrical tangent position. If the tangent
      // point is inside the covered latitude range, calculate also the 
      // geometrical altitude of the tangent point and the top of atmosphere.
      //
//      Array<double>  geom_tan_pos(0);
//      double geom_tan_z=-2, geom_tan_atmtop=-1;  // OK values if the variables

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = rte_pos[0];
      ppath.pos(0,1) = rte_pos[1];
      ppath.pos(0,2) = rte_pos[2];
      ppath.los(0,0) = rte_los[0];
      ppath.los(0,1) = rte_los[1];


     // Geometrical altitude
      ppath.z[0] = ppath.pos(0,0) - rv_geoid;

     // Use below the values in ppath (instead of rte_pos and rte_los) as 
     // they can be modified on the way.
     
     // Grid positions
      ppath.gp_lat[0].idx   = gp_lat.idx;
      ppath.gp_lat[0].fd[0] = gp_lat.fd[0];
      ppath.gp_lat[0].fd[1] = gp_lat.fd[1];
      ppath.gp_lon[0].idx   = gp_lon.idx;
      ppath.gp_lon[0].fd[0] = gp_lon.fd[0];
      ppath.gp_lon[0].fd[1] = gp_lon.fd[1];

      // Create a vector with the geometrical altitude of the pressure 
      // surfaces for the sensor latitude and use it to get ppath.gp_p.
      Vector z_grid(np);
      if ( z_field_is_1D )
        { z_grid=z_field(joker,0,0);}
      else
        { z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, 
                       gp_lat, gp_lon );}
      gridpos( ppath.gp_p, z_grid, ppath.z );


      // Handle possible numerical problems for grid positions
      gridpos_check_fd( ppath.gp_p[0] );
      gridpos_check_fd( ppath.gp_lat[0] );
      gridpos_check_fd( ppath.gp_lon[0] );

    }  // End 3D
}


//! cum_l_stepCalc

/*!
   Returns a vector of cumulative pathlengths for a given ppath by adding ppath.l_step

   \param cum_l_step    Output: vector of cumulative pathlengths.
   \param ppath         a Ppath.

   \author Cory Davis
   \date   2003-06-19
*/

void cum_l_stepCalc(
                      Vector& cum_l_step,
                      const Ppath& ppath
                      )
  {
    cum_l_step.resize(ppath.np);
    Numeric cumsum = 0.0;
    cum_l_step[0] = 0.0;
    for (Index i=0; i<ppath.np-1; i++)
      {
        cumsum += ppath.l_step[i];
        cum_l_step[i+1] = cumsum;
      }
  }             

//! findZ11max
/*! 
  The direction sampling method requires a bounding value for Z11.
  This returns a vector with the maximum value of Z11 for each particle type.

  \param[out] Z11maxvector Maximum value of Z11 for each particle type
  \param[in]  scat_data_mono

  \author Cory Davis
  \date   2004-31-1
*/
void findZ11max(Vector& Z11maxvector,
           const ArrayOfSingleScatteringData& scat_data_mono)
{
  Index np=scat_data_mono.nelem();
  Z11maxvector.resize(np);

  for(Index i = 0;i<np;i++)
    {
      switch(scat_data_mono[i].ptype){
      case PTYPE_MACROS_ISO:
        {
          Z11maxvector[i]=max(scat_data_mono[i].pha_mat_data(0,joker,joker,0,0,0,0));
        }
      case PTYPE_HORIZ_AL:
        {
          Z11maxvector[i]=max(scat_data_mono[i].pha_mat_data(0,joker,joker,0,joker,0,0));
        }
      default:
        Z11maxvector[i]=max(scat_data_mono[i].pha_mat_data(0,joker,joker,joker,joker,joker,0));
      }
    }
}




//! is_anyptype30
/*!
Some operations in Monte Carlo simulations are different depending on the 
particle type of the scattering particles.  This function searches 
scat_data_mono to determine if any of the particle types have ptype=30

\author Cory Davis
\date 2004-1-31

*/

bool is_anyptype30(const ArrayOfSingleScatteringData& scat_data_mono)
{
  Index np=scat_data_mono.nelem();
  bool anyptype30=false;
  Index i=0;
  while(i < np && anyptype30==false)
    {
      if(scat_data_mono[i].ptype==PTYPE_HORIZ_AL)
        {
          anyptype30=true;
        }
      i+=1;
    }
  return anyptype30;
}

//! iwp_cloud_opt_pathCalc
/*!

Calculated the ice water path(iwp) and cloud optical path (cloud_opt_path) for 
a given sensor_pos_ and sensor_los


\author Cory Davis
\date 2006-6-15

*/


void iwp_cloud_opt_pathCalc(
                            Numeric& iwp,
                            Numeric& cloud_opt_path,
                            //input
                            const Vector&         rte_pos,
                            const Vector&         rte_los,
                            const Agenda&         ppath_step_agenda,
                            const Vector&         p_grid,
                            const Vector&         lat_grid, 
                            const Vector&         lon_grid, 
                            const Matrix&         r_geoid, 
                            const Matrix&         z_surface,
                            const Tensor3&        z_field, 
                            const Tensor3&        t_field, 
                            const Tensor4&        vmr_field, 
                            const ArrayOfIndex&   cloudbox_limits, 
                            const Tensor4&        pnd_field,
                            const ArrayOfSingleScatteringData& scat_data_mono,
                            const Vector&          particle_masses
                            )
{
  //internal declarations
  Ppath ppath,ppath_step;
  Vector local_rte_pos=rte_pos;
  Vector local_rte_los=rte_los;
  iwp=0;
  cloud_opt_path=0;
  //calculate ppath to cloudbox boundary
  ppath_calc( ppath, ppath_step, ppath_step_agenda, 3, 
              p_grid, lat_grid, lon_grid, z_field, r_geoid, z_surface, 
              1, cloudbox_limits, local_rte_pos, local_rte_los, 1 );
  //if this ppath hit a cloud, now take ppath inside cloud
  if (ppath_what_background(ppath)>2)
    {
      //move to the intersection of the cloudbox boundary and the 
      //propagation path
      local_rte_pos=ppath.pos(ppath.np-1,joker);
      local_rte_los=ppath.los(ppath.np-1,joker);

      Range p_range(cloudbox_limits[0], 
                    cloudbox_limits[1]-cloudbox_limits[0]+1);
      Range lat_range(cloudbox_limits[2], 
                      cloudbox_limits[3]-cloudbox_limits[2]+1);
      Range lon_range(cloudbox_limits[4], 
                      cloudbox_limits[5]-cloudbox_limits[4]+1);

      ppath_calc( ppath, ppath_step, ppath_step_agenda, 3, 
                  p_grid, lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                  1, cloudbox_limits, local_rte_pos, local_rte_los, 0 );

      Matrix  pnd_ppath(particle_masses.nelem(),ppath.np);
      Vector t_ppath(ppath.np);
      Vector   p_ppath(ppath.np);
      Matrix   vmr_ppath(vmr_field.nbooks(),ppath.np);
      Matrix ext_mat_part(1, 1, 0.0);
      Vector abs_vec_part(1, 0.0);

      cloud_atm_vars_by_gp(p_ppath,t_ppath,vmr_ppath,pnd_ppath,ppath.gp_p,
                           ppath.gp_lat,ppath.gp_lon,cloudbox_limits,p_grid[p_range],
                           lat_grid[lat_range],lon_grid[lon_range],
                           t_field(p_range,lat_range,lon_range),
                           vmr_field(joker,p_range,lat_range,lon_range),
                           pnd_field);

      Vector k_vec(ppath.np);
      Vector iwc_vec(ppath.np);
      //calculate integrands
      for (Index i = 0; i < ppath.np ; ++i)
        {
          opt_propCalc(ext_mat_part,abs_vec_part,ppath.los(i,0),ppath.los(i,1),scat_data_mono,
                       1, pnd_ppath(joker, i), t_ppath[i]);
          k_vec[i]=ext_mat_part(0,0);
          Vector pnd_vec=pnd_ppath(joker, i);
          assert(pnd_vec.nelem()==particle_masses.nelem());
          iwc_vec[i]=pnd_vec*particle_masses;//hopefully this is the dot product
        }
      //integrate IWP and optical properties
      for (Index i = 0; i < ppath.np-1 ; ++i)
        {
          //Add to path integrals
          iwp+=0.5*(iwc_vec[i+1]+iwc_vec[i])*ppath.l_step[i];
          cloud_opt_path+=0.5*(k_vec[i+1]+k_vec[i])*ppath.l_step[i];
        }
    }
}


//! matrix_exp_p30
/*!
When we have p30 particles, and therefore the extinction matrix has a 
block diagonal form with 3 independent elements, the matrix expontential
can be calculated very quickly and exactly using this function.

\param M output matrix
\param A input matrix (must be of the form described above)

\author Cory Davis
\date 2005-3-2
*/
void matrix_exp_p30(MatrixView& M,
                    ConstMatrixView& A)
{
  Index m=A.nrows();
  assert( A.ncols()==m );
  M=0;
  Numeric a=A(0,0);
  Numeric b=A(0,1);
  M(0,0)=cosh(b);
  M(1,1)=cosh(b);
  M(0,1)=sinh(b);
  M(1,0)=sinh(b);
  if ( m>2 )
    {
      Numeric c=A(2,3);
      M(2,2)=cos(c);
      if ( m > 3 )
        {
          M(2,3)=sin(c);
          M(3,2)=-sin(c);
        }
    }
  M*=exp(a);    
}



//! mcPathTraceGeneral
/*!
Performs the tasks of pathlength sampling,
ray tracing (but now only as far as determined by pathlength
sampling) and calculation of the evolution operator and several 
atmospheric variables at the new point.

\author Cory Davis
\date 2005-2-21

*/

void mcPathTraceGeneral(MatrixView&           evol_op,
                        Vector&               abs_vec_mono,
                        Numeric&              temperature,
                        MatrixView&           ext_mat_mono,
                        Rng&                  rng,
                        Vector&               rte_pos,
                        Vector&               rte_los,
                        Vector&               pnd_vec,
                        Numeric&              g,
                        Ppath&                ppath_step,
                        Index&                termination_flag,
                        bool&                 inside_cloud,
                        //Numeric&              rte_pressure,
                        //Vector&               rte_vmr_list,
                        const Agenda&         opt_prop_gas_agenda,
                        const Agenda&         abs_scalar_gas_agenda,
                        const Index&          stokes_dim,
                        const Index&          f_index,
                        const Vector&         p_grid,
                        const Vector&         lat_grid,
                        const Vector&         lon_grid,
                        const Tensor3&        z_field,
                        const Matrix&         r_geoid,
                        const Matrix&         z_surface,
                        const Tensor3&        t_field,
                        const Tensor4&        vmr_field,
                        const ArrayOfIndex&   cloudbox_limits,
                        const Tensor4&        pnd_field,
                        const ArrayOfSingleScatteringData& scat_data_mono,
                        const Index&          z_field_is_1D)

{ 
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2),pnd_vecArray(2);
  Vector tArray(2);
  Vector cum_l_step;
  Matrix T(stokes_dim,stokes_dim);
  Numeric k;
  Numeric ds;
  Index np=0;
  Index   istep = 0;            // Counter for number of steps
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim,0.0);
  Matrix old_evol_op(stokes_dim,stokes_dim);
  //g=0;k=0;ds=0;
  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1]=evol_op;
  //initialise Ppath with ppath_start_stepping
  Index rubbish=z_field_is_1D;rubbish+=1;
  ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, r_geoid, z_surface,
                        0, cloudbox_limits, false, 
                        rte_pos, rte_los );
  //Use cloudbox_ppath_start_stepping to avoid unnecessary z_at_latlon calls.
  //cloudbox_ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
  //                               lon_grid, z_field, r_geoid, z_surface, rte_pos,
  //                               rte_los ,z_field_is_1D);
  Range p_range(cloudbox_limits[0], 
                cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range(cloudbox_limits[2], 
                cloudbox_limits[3]-cloudbox_limits[2]+1);
  Range lon_range(cloudbox_limits[4], 
                cloudbox_limits[5]-cloudbox_limits[4]+1);
  termination_flag=0;

  inside_cloud=is_inside_cloudbox( ppath_step, cloudbox_limits,true );
  
  if (inside_cloud)
    {
      cloudy_rt_vars_at_gp(ext_mat_mono,abs_vec_mono,pnd_vec,temperature,
                  opt_prop_gas_agenda,abs_scalar_gas_agenda,
                  stokes_dim, f_index, ppath_step.gp_p[0], ppath_step.gp_lat[0],
                  ppath_step.gp_lon[0],p_grid[p_range],lat_grid[lat_range], 
                  lon_grid[lon_range],t_field(p_range,lat_range,lon_range), 
                  vmr_field(joker,p_range,lat_range,lon_range),pnd_field,
                  scat_data_mono, cloudbox_limits,ppath_step.los(0,joker));
    }
  else
    {
      clear_rt_vars_at_gp( ext_mat_mono, abs_vec_mono, temperature, 
            opt_prop_gas_agenda, abs_scalar_gas_agenda, f_index, 
            ppath_step.gp_p[0], ppath_step.gp_lat[0], ppath_step.gp_lon[0],
            p_grid, lat_grid, lon_grid, t_field, vmr_field );
      pnd_vec=0.0;
    }
  ext_matArray[1]=ext_mat_mono;
  abs_vecArray[1]=abs_vec_mono;
  tArray[1]=temperature;
  pnd_vecArray[1]=pnd_vec;
  //draw random number to determine end point
  Numeric r = rng.draw();
  while ((evol_op(0,0)>r)&(!termination_flag))
    {

      istep++;
      //we consider more than 5000
      // path points to be an indication on that the calcululations have
      // got stuck in an infinite loop.
      if( istep > 5000 )
        {
          throw runtime_error(
                            "5000 path points have been reached. Is this an infinite loop?" );
        }
      evol_opArray[0]=evol_opArray[1];
      ext_matArray[0]=ext_matArray[1];
      abs_vecArray[0]=abs_vecArray[1];
      tArray[0]=tArray[1];
      pnd_vecArray[0]=pnd_vecArray[1];
      //perform single path step using ppath_step_geom_3d
      ppath_step_geom_3d(ppath_step, p_grid, lat_grid, lon_grid, z_field,
                         r_geoid, z_surface, -1);
      np=ppath_step.np;//I think this should always be 2.
      cum_l_stepCalc(cum_l_step, ppath_step);
      //path_step should now have two elements.
      //calculate evol_op
      inside_cloud=is_inside_cloudbox( ppath_step, cloudbox_limits, true );
      if (inside_cloud)
        {
          cloudy_rt_vars_at_gp(ext_mat_mono,abs_vec_mono,pnd_vec,temperature,
               opt_prop_gas_agenda,abs_scalar_gas_agenda, stokes_dim, f_index,
               ppath_step.gp_p[np-1],ppath_step.gp_lat[np-1],
               ppath_step.gp_lon[np-1],p_grid[p_range], lat_grid[lat_range], 
               lon_grid[lon_range],t_field(p_range,lat_range,lon_range), 
               vmr_field(joker,p_range,lat_range,lon_range),pnd_field,
               scat_data_mono, cloudbox_limits,ppath_step.los(np-1,joker));
        }
      else
        {
          clear_rt_vars_at_gp(ext_mat_mono,abs_vec_mono,temperature, 
               opt_prop_gas_agenda, abs_scalar_gas_agenda, f_index, 
               ppath_step.gp_p[np-1], ppath_step.gp_lat[np-1],
               ppath_step.gp_lon[np-1], p_grid, lat_grid, lon_grid, t_field, 
                                                                    vmr_field);
          pnd_vec=0.0;
        }
      ext_matArray[1]=ext_mat_mono;
      abs_vecArray[1]=abs_vec_mono;
      tArray[1]=temperature;
      pnd_vecArray[1]=pnd_vec;
      opt_depth_mat=ext_matArray[1];
      opt_depth_mat+=ext_matArray[0];
      opt_depth_mat*=-cum_l_step[np-1]/2;
      incT=0;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index j=0;j<stokes_dim;j++)
            {
              incT(j,j)=exp(opt_depth_mat(j,j));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[0],incT);
      evol_opArray[1]=evol_op;
     
      //check whether hit ground or space
      Index bg = ppath_what_background(ppath_step);
      if ( bg==2 )
        {
          //we have hit the surface
          termination_flag=2;
          
        }
      //else if ( is_gridpos_at_index_i( ppath_step.gp_p[np-1], p_grid.nelem() - 1 ) )
      else if (fractional_gp(ppath_step.gp_p[np-1])>=(p_grid.nelem() - 1)-1e-3)
        {
          //we have left the top of the atmosphere
          termination_flag=1;
          
        }


    }
  if (termination_flag)
    {//we must have reached the cloudbox boundary
      //
      rte_pos = ppath_step.pos(np-1,Range(0,3));
      rte_los = ppath_step.los(np-1,joker);
      g=evol_op(0,0);
    }
  else
    {
      //find position...and evol_op..and everything else required at the new
      //scattering/emission point
      k=-log(incT(0,0))/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
      ds=log(evol_opArray[0](0,0)/r)/k;
      g=k*r;
      Vector x(2,0.0);
      //interpolate atmospheric variables required later
      //length 2
      ArrayOfGridPos gp(1);
      x[1]=cum_l_step[np-1];
      Vector itw(2);
  
      gridpos(gp, x, ds);
      assert(gp[0].idx==0);
      interpweights(itw,gp[0]);
      ext_mat_mono = interp(itw,ext_matArray,gp[0]);
      opt_depth_mat = ext_mat_mono;
      opt_depth_mat+=ext_matArray[gp[0].idx];
      opt_depth_mat*=-ds/2;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index i=0;i<stokes_dim;i++)
            {
              incT(i,i)=exp(opt_depth_mat(i,i));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[gp[0].idx],incT);
      abs_vec_mono = interp(itw, abs_vecArray,gp[0]);
      temperature=interp(itw,tArray,gp[0]);
      pnd_vec=interp(itw,pnd_vecArray,gp[0]);
      if (np > 2)
        {
          gridpos(gp,cum_l_step,ds);
          interpweights(itw,gp[0]);
        }
      for (Index i=0; i<2; i++)
        {
          rte_pos[i] = interp(itw,ppath_step.pos(Range(joker),i),gp[0]);
          rte_los[i] = interp(itw,ppath_step.los(Range(joker),i),gp[0]);
        }
      rte_pos[2] = interp(itw,ppath_step.pos(Range(joker),2),gp[0]);
    }
}

//! mcPathTraceIPA
/*!
Performs the tasks of pathlength sampling,
ray tracing (but now only as far as determined by pathlength
sampling) and calculation of the evolution operator and several 
atmospheric variables at the new point.  THis is the same as mcPathTraceGeneral
modified for the independent pixel approximation.

\author Cory Davis
\date 2005-2-21

*/

void mcPathTraceIPA(MatrixView&           evol_op,
                    Vector&               abs_vec_mono,
                    Numeric&              temperature,
                    MatrixView&           ext_mat_mono,
                    Rng&                  rng,
                    Vector&               rte_pos,
                    Vector&               rte_los,
                    Vector&               pnd_vec,
                    Numeric&              g,
                    Ppath&                ppath_step,
                    Index&                termination_flag,
                    bool&                 inside_cloud,
                    const Agenda&         opt_prop_gas_agenda,
                    const Agenda&         abs_scalar_gas_agenda,
                    const Index&          stokes_dim,
                    const Index&          f_index,
                    const Vector&         p_grid,
                    const Vector&         lat_grid,
                    const Vector&         lon_grid,
                    const Tensor3&        z_field,
                    const Matrix&         r_geoid,
                    const Matrix&         z_surface,
                    const Tensor3&        t_field,
                    const Tensor4&        vmr_field,
                    const ArrayOfIndex&   cloudbox_limits,
                    const Tensor4&        pnd_field,
                    const ArrayOfSingleScatteringData& scat_data_mono,
                    const Index&          z_field_is_1D,
                    const Ppath&          ppath)

{ 

  //Internal declarations
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2),pnd_vecArray(2);
  GridPos gp_p,gp_lat,gp_lon;
  Vector    itw(4);
  Vector tArray(2);
  Vector cum_l_step;
  Vector z_grid(p_grid.nelem());
  Matrix T(stokes_dim,stokes_dim);
  Numeric k;
  Numeric x,y,z;
  Numeric ds,lstep=0.,dx,dy,dz;
  Index   istep = 0;            // Counter for number of steps
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim,0.0);
  Matrix old_evol_op(stokes_dim,stokes_dim);
  Numeric rv_geoid=0.;
  Numeric rv_surface=0.;
  Numeric alt;
  Numeric gpnum;
  const Index np_ipa=ppath.np;
  Index i_closest=0;
  Numeric gp_diff;
  Vector lon_iw(2),lat_iw(2);
  Vector l_vec(2); 
  GridPos gp;
  Vector itw1d(2);
        
  //Initialisations
  
  if (z_field_is_1D)
    {
      z_grid=z_field(joker,0,0);
      rv_geoid=r_geoid(0,0);
      rv_surface=rv_geoid+z_surface(0,0);
    }

  id_mat(evol_op);
  evol_opArray[1]=evol_op;
  //initialise Ppath with ppath_start_stepping
  ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, r_geoid, z_surface,
                        0, cloudbox_limits, false, 
                        rte_pos, rte_los );

  gp_p=ppath_step.gp_p[0];
  gp_lat=ppath_step.gp_lat[0];
  gp_lon=ppath_step.gp_lon[0];
  rte_pos=ppath_step.pos(0,joker);
  rte_los=ppath_step.los(0,joker);

  Range p_range(cloudbox_limits[0], 
                cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range(cloudbox_limits[2], 
                cloudbox_limits[3]-cloudbox_limits[2]+1);
  Range lon_range(cloudbox_limits[4], 
                cloudbox_limits[5]-cloudbox_limits[4]+1);
  termination_flag=0;

  inside_cloud=is_inside_cloudbox( ppath_step, cloudbox_limits,false );
  
  if (inside_cloud)
    {
      cloudy_rt_vars_at_gp(ext_mat_mono,abs_vec_mono,pnd_vec,temperature,
              opt_prop_gas_agenda,abs_scalar_gas_agenda, stokes_dim, f_index,
              gp_p, gp_lat, gp_lon,p_grid[p_range], 
              lat_grid[lat_range], lon_grid[lon_range], 
              t_field(p_range,lat_range,lon_range), 
              vmr_field(joker,p_range,lat_range,lon_range),
              pnd_field,scat_data_mono, cloudbox_limits,rte_los );
    }
  else
    {
      clear_rt_vars_at_gp( ext_mat_mono,abs_vec_mono,temperature, 
                           opt_prop_gas_agenda, abs_scalar_gas_agenda, f_index,
                           gp_p, gp_lat, gp_lon,
                           p_grid, lat_grid, lon_grid, t_field, vmr_field);
      pnd_vec=0.0;
    }
  ext_matArray[1]=ext_mat_mono;
  abs_vecArray[1]=abs_vec_mono;
  tArray[1]=temperature;
  pnd_vecArray[1]=pnd_vec;
  //draw random number to determine end point
  Numeric r = rng.draw();
  while ((evol_op(0,0)>r)&(!termination_flag))
    {
      //check we are not in an infinite loop (assert for ease of debugging
      istep++;
      assert(istep<=5000);
      
      //these array variables hold values from the two ends of a path step.
      //Here we make the values for the last point of the previous step, the values
      //for the first poinst of the current step.
      evol_opArray[0]=evol_opArray[1];
      ext_matArray[0]=ext_matArray[1];
      abs_vecArray[0]=abs_vecArray[1];
      tArray[0]=tArray[1];
      pnd_vecArray[0]=pnd_vecArray[1];
      
      //Choose sensible path step length.  This is done by considering the dimensions of the
      //current grid cell and the LOS direction. The aim is to move only one grid cell.
      //dy=cos(DEG2RAD*rte_los[1])*sin(DEG2RAD*rte_los[0])*
      //  (lat_grid[gp_lat.idx+1]-lat_grid[gp_lat.idx])*DEG2RAD*r_geoid(gp_lat.idx,gp_lon.idx);
      //dx=sin(DEG2RAD*rte_los[1])*sin(DEG2RAD*rte_los[0])*
      //  (lon_grid[gp_lon.idx+1]-lon_grid[gp_lon.idx])*DEG2RAD*r_geoid(gp_lat.idx,gp_lon.idx)*
      //  cos(DEG2RAD*rte_pos[1]);
      //dz=cos(DEG2RAD*rte_los[0])*(z_field(gp_p.idx+1,gp_lat.idx,gp_lon.idx)-
      //                            z_field(gp_p.idx,gp_lat.idx,gp_lon.idx));
      //lstep=sqrt(dx*dx+dy*dy+dz*dz);
      //take the minimum grid dimension as the path step length
      Numeric Dy=(lat_grid[gp_lat.idx+1]-lat_grid[gp_lat.idx])*DEG2RAD*
        r_geoid(gp_lat.idx,gp_lon.idx);
      Numeric Dx=(lon_grid[gp_lon.idx+1]-lon_grid[gp_lon.idx])*
        DEG2RAD*r_geoid(gp_lat.idx,gp_lon.idx)*cos(DEG2RAD*rte_pos[1]);
      Numeric Dz=z_field(gp_p.idx+1,gp_lat.idx,gp_lon.idx)-
        z_field(gp_p.idx,gp_lat.idx,gp_lon.idx);
      Numeric Dxy=sqrt(Dx*Dx+Dy*Dy);
      if (Dxy/abs(tan(rte_los[0]*DEG2RAD))>Dz){dz=Dz;}
      else {dz=Dxy/abs(tan(rte_los[0]*DEG2RAD));}
      Numeric dxy;
      if(dz*abs(tan(rte_los[0]*DEG2RAD))>Dxy){dxy=Dxy;}
      else{dxy=dz*abs(tan(rte_los[0]*DEG2RAD));}
      lstep=sqrt(dxy*dxy+dz*dz);
      //calculate new position
      //current position and direction vector
      poslos2cart( x, y, z, dx, dy, dz, rte_pos[0], rte_pos[1], rte_pos[2], 
                   rte_los[0], rte_los[1] );
      //new_position
      x+=lstep*dx;
      y+=lstep*dy;
      z+=lstep*dz;
      //back to spherical coords
      cart2poslos(rte_pos[0],rte_pos[1],rte_pos[2],rte_los[0],rte_los[1],
                  x,y,z,dx,dy,dz);
      //get new grid_positions
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      gridpos( gp_lon, lon_grid, rte_pos[2] );
      if (!z_field_is_1D)
        {
          z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, gp_lat, gp_lon );
          interpweights( itw, gp_lat, gp_lon );
          rv_geoid  = interp( itw, r_geoid, gp_lat, gp_lon );
          rv_surface = rv_geoid + interp( itw, z_surface, gp_lat, gp_lon );
        }
      alt = rte_pos[0] - rv_geoid;
      //Have we left the top of the atmosphere?
      if ((alt>=z_grid[p_grid.nelem()-1])&(rte_los[0]<=90))
        {
          termination_flag=1;
          //shift relavent position variables to the appropriate point
          //z is linearly interpolated with respect to lat and lon
          //I don't think this is overly important.  Just change gp_p and rte_pos
          alt=z_grid[p_grid.nelem()-1];
          rte_pos[0]=alt+rv_geoid;
        }
      //Have we hit the surface?
      if( (rte_pos[0] <= rv_surface)&(rte_los[0]>=90) )
        {
          termination_flag=2;
          rte_pos[0]=rv_surface;
          alt = rte_pos[0] - rv_geoid;
        }
      gridpos( gp_p, z_grid, alt );

      //now project the new point back to the original IPA path
      //find lat and lon gridpoints on reference ppath_ipa 
      gpnum=fractional_gp(gp_p);
      //search through ppath_ipa to find the closest grid point
      gp_diff=abs(fractional_gp(ppath.gp_p[0])-gpnum);
      for (Index i=1;i<np_ipa;i++)
        {
          Numeric new_gp_diff=abs(fractional_gp(ppath.gp_p[i])-gpnum);
          if (new_gp_diff<=gp_diff)
            {
              gp_diff=new_gp_diff;
              i_closest=i;
            }
          else
            {
              //getting further away - can stop now
              break;
            }
        }
      gp_lat=ppath.gp_lat[i_closest];
      gp_lon=ppath.gp_lon[i_closest];
      interpweights(lon_iw,gp_lon);
      interpweights(lat_iw,gp_lat);
      //      
      if (!z_field_is_1D)
       {
          interpweights( itw, gp_lat, gp_lon );
          rv_geoid  = interp( itw, r_geoid, gp_lat, gp_lon );          
        }
      rte_pos[0]=rv_geoid+interp_atmfield_by_gp(3,p_grid,lat_grid,lon_grid,z_field,
                                                "z_field",gp_p,gp_lat,gp_lon);
      rte_pos[1]=interp(lat_iw, lat_grid,gp_lat);
      rte_pos[2]=interp(lon_iw, lon_grid,gp_lon);

      //calculate evolution operator
      inside_cloud=is_gp_inside_cloudbox(gp_p,gp_lat,gp_lon,
                                         cloudbox_limits, false );
      //calculate RT variables at new point
      if (inside_cloud)
        {
          cloudy_rt_vars_at_gp(ext_mat_mono, abs_vec_mono, pnd_vec, temperature,
                   opt_prop_gas_agenda, abs_scalar_gas_agenda, stokes_dim, 
                   f_index, gp_p,gp_lat, gp_lon, p_grid[p_range], 
                   lat_grid[lat_range], lon_grid[lon_range], 
                   t_field(p_range,lat_range,lon_range), 
                   vmr_field(joker,p_range,lat_range,lon_range),
                   pnd_field, scat_data_mono, cloudbox_limits,rte_los);
        }
      else
        {
          clear_rt_vars_at_gp( ext_mat_mono, abs_vec_mono, temperature, 
                      opt_prop_gas_agenda, abs_scalar_gas_agenda, f_index,
                      gp_p, gp_lat, gp_lon, p_grid, lat_grid, lon_grid, t_field, 
                                                                     vmr_field);
          pnd_vec=0.0;
        }
      //put these variables in the last Array slot
      ext_matArray[1]=ext_mat_mono;
      abs_vecArray[1]=abs_vec_mono;
      tArray[1]=temperature;
      pnd_vecArray[1]=pnd_vec;
      opt_depth_mat=ext_matArray[1];
      //now calculate evolution operator
      opt_depth_mat+=ext_matArray[0];
      opt_depth_mat*=-lstep/2;
      incT=0;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index j=0;j<stokes_dim;j++)
            {
              incT(j,j)=exp(opt_depth_mat(j,j));
            }
        }
      else
        {
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[0],incT);
      assert(incT(0,0)<1);
      assert(evol_op(0,0)<1);
      evol_opArray[1]=evol_op;

    }
  if (termination_flag)
    {
      //we must have reached the surface or the TOA
      g=evol_op(0,0);
    }
  else
    {
      //then we have met the monte carlo pathlength criteria somewhere inside 
      //the atmosphere and between the two points corresponding to the *Array 
      //variables.

      //find position corresponding to the pathlength criteria
      k=-log(incT(0,0))/lstep;
      //length of path segment required by critera
      ds=log(evol_opArray[0](0,0)/r)/k;
      g=k*r;
      l_vec=0.0;
      l_vec[1]=lstep;
      //gridpos and interpolations for required step length
      gridpos(gp, l_vec, ds);
      interpweights(itw1d,gp);
      //interpolate optical properties at required point
      ext_mat_mono = interp(itw1d,ext_matArray,gp);
      opt_depth_mat = ext_mat_mono;
      opt_depth_mat+=ext_matArray[gp.idx];
      opt_depth_mat*=-ds/2;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index i=0;i<stokes_dim;i++)
            {
              incT(i,i)=exp(opt_depth_mat(i,i));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[gp.idx],incT);
      abs_vec_mono = interp(itw1d, abs_vecArray,gp);
      temperature=interp(itw1d,tArray,gp);
      pnd_vec=interp(itw1d,pnd_vecArray,gp);
      //move cartesion coordinates back to required point.
      x+=(ds-lstep)*dx;
      y+=(ds-lstep)*dy;
      z+=(ds-lstep)*dz;
      //and convert back to spherical.
      cart2poslos(rte_pos[0],rte_pos[1],rte_pos[2],rte_los[0],rte_los[1],x,y,z,dx,dy,dz);
      
    }
}

  

//! mcPathTrace
/*!
Performs the tasks of pathlength sampling,
ray tracing (but now only as far as determined by pathlength
sampling) and calculation of the evolution operator and several 
atmospheric variables at the new point.

\author Cory Davis
\date 2005-2-21

*/

void mcPathTrace(MatrixView&           evol_op,
                 VectorView&           abs_vec_mono,
                 Numeric&              rte_temperature,
                 MatrixView&           ext_mat_mono,
                 Rng&                  rng,
                 Vector&               rte_pos,
                 Vector&               rte_los,
                 Vector&               pnd_vec,
                 Numeric&              g,
                 bool&                 left_cloudbox,
                 const Agenda&         opt_prop_gas_agenda,
                 const Agenda&         abs_scalar_gas_agenda,
                 const Index&          stokes_dim,
                 const Index&          f_index,
                 const Vector&         p_grid,
                 const Vector&         lat_grid,
                 const Vector&         lon_grid,
                 const Tensor3&        z_field,
                 const Matrix&         r_geoid,
                 const Matrix&         z_surface,
                 const Tensor3&        t_field,
                 const Tensor4&        vmr_field,
                 const ArrayOfIndex&   cloudbox_limits,
                 const Tensor4&        pnd_field,
                 const ArrayOfSingleScatteringData& scat_data_mono,
                 const Index&          z_field_is_1D)

{
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2),pnd_vecArray(2);
  Vector tArray(2);
  Vector cum_l_step;
  Matrix T(stokes_dim,stokes_dim);
  Ppath ppath_step;
  Numeric k;
  Numeric ds;
  Index np=0;
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim,0.0);
  Matrix old_evol_op(stokes_dim,stokes_dim);
  left_cloudbox=false;
  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1]=evol_op;
  //initialise Ppath with ppath_start_stepping
  /*ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, r_geoid, z_surface,
                        1, cloudbox_limits, false, 
                        rte_pos, rte_los );*/
  //Use cloudbox_ppath_start_stepping to avoid unnecessary z_at_latlon calls.
  cloudbox_ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                                 lon_grid, z_field, r_geoid, z_surface, rte_pos,
                                 rte_los ,z_field_is_1D);
  Range p_range(cloudbox_limits[0], 
                cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range(cloudbox_limits[2], 
                cloudbox_limits[3]-cloudbox_limits[2]+1);
  Range lon_range(cloudbox_limits[4], 
                cloudbox_limits[5]-cloudbox_limits[4]+1);
  cloudy_rt_vars_at_gp( ext_mat_mono, abs_vec_mono,pnd_vec, rte_temperature,
              opt_prop_gas_agenda, abs_scalar_gas_agenda, stokes_dim, f_index, 
              ppath_step.gp_p[0], ppath_step.gp_lat[0], ppath_step.gp_lon[0],
              p_grid[p_range], lat_grid[lat_range], lon_grid[lon_range], 
              t_field(p_range,lat_range,lon_range),
              vmr_field(joker,p_range,lat_range,lon_range), pnd_field,
              scat_data_mono, cloudbox_limits, ppath_step.los(0,joker) );
  ext_matArray[1]=ext_mat_mono;
  abs_vecArray[1]=abs_vec_mono;
  tArray[1]=rte_temperature;
  pnd_vecArray[1]=pnd_vec;
  //draw random number to determine end point
  Numeric r = rng.draw();
  while ((evol_op(0,0)>r)&(!left_cloudbox))
    {
      evol_opArray[0]=evol_opArray[1];
      ext_matArray[0]=ext_matArray[1];
      abs_vecArray[0]=abs_vecArray[1];
      tArray[0]=tArray[1];
      pnd_vecArray[0]=pnd_vecArray[1];
      //perform single path step using ppath_step_geom_3d
      ppath_step_geom_3d(ppath_step, p_grid, lat_grid, lon_grid, z_field,
                         r_geoid, z_surface, -1);
      np=ppath_step.np;//I think this should always be 2.
      cum_l_stepCalc(cum_l_step, ppath_step);
      //path_step should now have two elements.
      //calculate evol_op
      cloudy_rt_vars_at_gp( ext_mat_mono, abs_vec_mono, pnd_vec,
             rte_temperature, opt_prop_gas_agenda, abs_scalar_gas_agenda, 
             stokes_dim, f_index, ppath_step.gp_p[np-1], ppath_step.gp_lat[np-1],
             ppath_step.gp_lon[np-1], p_grid[p_range], lat_grid[lat_range], 
             lon_grid[lon_range], t_field(p_range,lat_range,lon_range), 
             vmr_field(joker,p_range,lat_range,lon_range),
             pnd_field,scat_data_mono, cloudbox_limits,
             ppath_step.los(np-1,joker));
      ext_matArray[1]=ext_mat_mono;
      abs_vecArray[1]=abs_vec_mono;
      tArray[1]=rte_temperature;
      pnd_vecArray[1]=pnd_vec;
      opt_depth_mat=ext_matArray[1];
      opt_depth_mat+=ext_matArray[0];
      opt_depth_mat*=-cum_l_step[np-1]/2;
      incT=0;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index j=0;j<stokes_dim;j++)
            {
              incT(j,j)=exp(opt_depth_mat(j,j));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[0],incT);
      evol_opArray[1]=evol_op;
      left_cloudbox=not(is_inside_cloudbox(ppath_step,cloudbox_limits,false));
      }

  if (left_cloudbox)
    {//we must have reached the cloudbox boundary
      //
      rte_pos = ppath_step.pos(np-1,Range(0,3));
      rte_los = ppath_step.los(np-1,joker);
      g=evol_op(0,0);
    }
  else
    {
      //find position...and evol_op..and everything else required at the new
      //scattering/emission point
      k=-log(incT(0,0))/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
      ds=log(evol_opArray[0](0,0)/r)/k;
      g=k*r;
      Vector x(2,0.0);
      //interpolate atmospheric variables required later
      //length 2
      ArrayOfGridPos gp(1);
      x[1]=cum_l_step[np-1];
      Vector itw(2);
  
      gridpos(gp, x, ds);
      assert(gp[0].idx==0);
      interpweights(itw,gp[0]);
      ext_mat_mono = interp(itw,ext_matArray,gp[0]);
      opt_depth_mat = ext_mat_mono;
      opt_depth_mat+=ext_matArray[gp[0].idx];
      opt_depth_mat*=-ds/2;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index i=0;i<stokes_dim;i++)
            {
              incT(i,i)=exp(opt_depth_mat(i,i));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[gp[0].idx],incT);
      abs_vec_mono = interp(itw, abs_vecArray,gp[0]);
      rte_temperature=interp(itw,tArray,gp[0]);
      pnd_vec=interp(itw,pnd_vecArray,gp[0]);
      if (np > 2)
        {
          gridpos(gp,cum_l_step,ds);
          interpweights(itw,gp[0]);
        }
      for (Index i=0; i<2; i++)
        {
          rte_pos[i] = interp(itw,ppath_step.pos(Range(joker),i),gp[0]);
          rte_los[i] = interp(itw,ppath_step.los(Range(joker),i),gp[0]);
        }
      rte_pos[2] = interp(itw,ppath_step.pos(Range(joker),2),gp[0]);
    }
}


//!  montecarloGetIncoming 

/*!
   Gets incoming radiance at the cloud box boundary in a single propagation 
direction, determined by the cloudbox Ppath 'pptahcloud'. 
Used in ScatteringMonteCarlo.  

   \author Cory Davis
   \date   2003-11-28
*/

void montecarloGetIncoming(
                           Matrix&               iy,
                           Vector&               rte_pos,
                           Vector&               rte_los,
                           Ppath&                ppath,
                           Ppath&                ppath_step,
                           const Agenda&         ppath_step_agenda,
                           const Agenda&         rte_agenda,
                           const Agenda&         iy_space_agenda,
                           const Agenda&         surface_prop_agenda,
                           const Agenda&         iy_cloudbox_agenda,
                           const Vector&         p_grid,
                           const Vector&         lat_grid,
                           const Vector&         lon_grid,
                           const Tensor3&        z_field,
                           const Tensor3&        t_field,
                           const Tensor4&        vmr_field,
                           const Matrix&         r_geoid,
                           const Matrix&         z_surface,
                           const ArrayOfIndex&   cloudbox_limits,
                           const Index&          atmosphere_dim,
                           const Vector&         f_grid,
                           const Index&          stokes_dim
                           )

{
  //call rte_calc without input checking, sensor stuff, or verbosity
  // Assign dummies for variables associated with sensor.
  Vector   mblock_za_grid_dummy(1);
  mblock_za_grid_dummy[0] = 0;
  Vector   mblock_aa_grid_dummy(0);
  //Index    antenna_dim_dummy = 1;
  Sparse   sensor_response_dummy;
  
  // Dummy for measurement vector
  Vector   y_dummy(0);
  
  const Index cloudbox_on_dummy=0;
  Tensor7 scat_i_p_dummy;
  Tensor7 scat_i_lat_dummy;
  Tensor7 scat_i_lon_dummy;
  
  Vector pos = rte_pos;
  Vector los = rte_los;
  iy_calc_no_jacobian( iy, ppath, ppath_step,
                       ppath_step_agenda, rte_agenda, iy_space_agenda,
                       surface_prop_agenda, iy_cloudbox_agenda,
                       atmosphere_dim, p_grid, lat_grid, lon_grid, z_field,
                       t_field, vmr_field, r_geoid, z_surface,
                       cloudbox_on_dummy, cloudbox_limits,
                       pos, los, f_grid,stokes_dim, 1);
  
  for (Index i = 0;i<stokes_dim;i++){assert(!isnan(iy(0,i)));}
}

//! opt_depth_calc  

/*!
   Used at the beginning of a MC simulation to calculate the optical depth
   between the sensor and the cloudbox boundary.

   \author Cory Davis
   \date   Late 2004
*/


Numeric opt_depth_calc(
                       Tensor3&   ext_mat,
                       Matrix&    abs_vec,
                       Numeric&   rte_pressure,
                       Numeric&   rte_temperature,
                       Vector&    rte_vmr_list,
                       const Ppath&     ppath,
                       const Agenda& opt_prop_gas_agenda,
                       const Agenda& abs_scalar_gas_agenda,
                       const Index&     f_index,
                       const Vector&    p_grid,
                       const Vector&    lat_grid,
                       const Vector&    lon_grid,
                       const Tensor3&   t_field,
                       const Tensor4&   vmr_field,
                       const Index&     atmosphere_dim)
{
  const Index np=ppath.np;  
  const Index   ns = vmr_field.nbooks();
  Vector t_ppath(np);
                
  Vector kvector(np);
  Numeric opt_depth=0;
  // Determine the pressure at each propagation path point
  Vector   p_ppath(np);
  Matrix   itw_p(np,2);
  //
  interpweights( itw_p, ppath.gp_p );      
  itw2p( p_ppath, p_grid, ppath.gp_p, itw_p );
  
  // Determine the atmospheric temperature and species VMR at 
  // each propagation path point
  Matrix   vmr_ppath(ns,np), itw_field;
  //
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid, lat_grid, 
                          lon_grid, ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
  //
  interp_atmfield_by_itw( t_ppath,  atmosphere_dim, p_grid, lat_grid, 
                          lon_grid, t_field, "t_field", 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
  // 
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( vmr_ppath(is, joker), atmosphere_dim,
                              p_grid, lat_grid, lon_grid, 
                              vmr_field(is, joker, joker,  joker), 
                              "vmr_field", ppath.gp_p, ppath.gp_lat, 
                              ppath.gp_lon, itw_field );
    }
  //Create array of extinction cefficients corresponding to each point in ppath
  for (Index i=0; i<ppath.np; i++)
    { 
      Matrix abs_scalar_gas;
      rte_pressure    = p_ppath[i];
      rte_temperature = t_ppath[i];
      rte_vmr_list    = vmr_ppath(joker,i);
      abs_scalar_gas_agendaExecute (abs_scalar_gas, f_index,
                                           rte_pressure, rte_temperature,
                                           rte_vmr_list,
                                           abs_scalar_gas_agenda, true);
      opt_prop_gas_agendaExecute (ext_mat, abs_vec, abs_scalar_gas,
                                  opt_prop_gas_agenda, true);
      kvector[i]=ext_mat(0,0,0);
    }
  for (Index i=1; i<ppath.np;i++)
    {
      opt_depth+=(kvector[i-1]+kvector[i])*ppath.l_step[i-1]/2;
    }
  return opt_depth;
}

//! opt_propCalc
/*!
Returns the extinction matrix and absorption vector due to scattering particles
from scat_data_mono

   \param ext_mat_mono               Output: extinction matrix
   \param abs_vec_mono           Output: absorption coefficient vector
   \param za              zenith angle of propagation direction
   \param aa              azimuthal angle of propagation
   \param scat_data_mono  workspace variable
   \param stokes_dim     workspace variable
   \param pnd_vec         vector pf particle number densities (one element per particle type)
   \param rte_temperature loacl temperature (workspace variable)

   \author Cory Davis
   \date   2004-7-16
*/

void opt_propCalc(
                  MatrixView& ext_mat_mono,
                  VectorView& abs_vec_mono,
                  const Numeric za,
                  const Numeric aa,
                  const ArrayOfSingleScatteringData& scat_data_mono,
                  const Index&          stokes_dim,
                  const VectorView& pnd_vec,
                  const Numeric& rte_temperature
                  )
{
  const Index N_pt = scat_data_mono.nelem();

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }
  Matrix ext_mat_mono_spt(stokes_dim,stokes_dim);
  Vector abs_vec_mono_spt(stokes_dim);

  ext_mat_mono=0.0;
  abs_vec_mono=0.0;  

  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      if (pnd_vec[i_pt]>0)
        {
          opt_propExtract(ext_mat_mono_spt,abs_vec_mono_spt,scat_data_mono[i_pt],za,aa,
                          rte_temperature,stokes_dim);
          ext_mat_mono_spt*=pnd_vec[i_pt];
          abs_vec_mono_spt*=pnd_vec[i_pt];
          ext_mat_mono+=ext_mat_mono_spt;
          abs_vec_mono+=abs_vec_mono_spt;
        }
    }

}

//! Extract ext_mat_mono and abs_vec_mono from a monochromatic SingleScatteringData object
/*!
  Given a monochromatic SingleScatteringData object, propagation directions, 
  and the temperature, this function returns the extinction matrix and the
  absorption coefficient vector
  Output:
  \param ext_mat_mono_spt the extinction matrix (spt: this is for a single particle)
  \param abs_vec_mono_spt the absorption coefficient vector
  Input:
  \param scat_data a monochromatic SingleScatteringData object
  \param za
  \param aa
  \param rte_temperature
  \param stokes_dim

  \author Cory Davis
  \date 2004-07-16

*/

void opt_propExtract(
                     MatrixView& ext_mat_mono_spt,
                     VectorView& abs_vec_mono_spt,
                     const SingleScatteringData& scat_data,
                     const Numeric& za,
                     const Numeric& aa,
                     const Numeric& rte_temperature,
                     const Index& stokes_dim
                     )
{

  GridPos t_gp;
  //interpolate over temperature
  gridpos(t_gp, scat_data.T_grid, rte_temperature);
  Numeric x=aa;//just to avoid warnings until we use ptype=10
  x+=1;        //
  switch (scat_data.ptype){

  case PTYPE_GENERAL:
    // This is only included to remove warnings about unused variables 
    // during compilation

    out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
    break;
    
  case PTYPE_MACROS_ISO:
    {
      assert (scat_data.ext_mat_data.ncols() == 1);
      
      Vector itw(2);
      interpweights(itw, t_gp);
      // In the case of randomly oriented particles the extinction matrix is 
      // diagonal. The value of each element of the diagonal is the
      // extinction cross section, which is stored in the database.
     
      ext_mat_mono_spt = 0.;
      
      ext_mat_mono_spt(0,0) = interp(itw,scat_data.ext_mat_data(0,joker,0,0,0),t_gp);
      
      abs_vec_mono_spt = 0;

      abs_vec_mono_spt[0] = interp(itw,scat_data.abs_vec_data(0,joker,0,0,0),t_gp);      
      
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

  case PTYPE_HORIZ_AL://Added by Cory Davis 9/12/03
    {
      assert (scat_data.ext_mat_data.ncols() == 3);
      
      // In the case of horizontally oriented particles the extinction matrix
      // has only 3 independent non-zero elements ext_mat_monojj, K12=K21, and K34=-K43.
      // These values are dependent on the zenith angle of propagation. The 
      // data storage format also makes use of the fact that in this case
      //K(za_sca)=K(180-za_sca). 

      // 1st interpolate data by za_sca
      GridPos za_gp;
      Vector itw(4);
      Numeric Kjj;
      Numeric K12;
      Numeric K34;
      
      if (za>90)
        {
          gridpos(za_gp,scat_data.za_grid,180-za);
        }
      else
        {
          gridpos(za_gp,scat_data.za_grid,za);
        }

      interpweights(itw,t_gp,za_gp);
      
      ext_mat_mono_spt=0.0;
      abs_vec_mono_spt=0.0;
      Kjj=interp(itw,scat_data.ext_mat_data(0,joker,joker,0,0),t_gp,za_gp);
      ext_mat_mono_spt(0,0)=Kjj;
      abs_vec_mono_spt[0] = interp(itw,scat_data.abs_vec_data(0,joker,joker,0,0),
                            t_gp,za_gp);
      if( stokes_dim == 1 ){
        break;
      }
      
      K12=interp(itw,scat_data.ext_mat_data(0,joker,joker,0,1),t_gp,za_gp);
      ext_mat_mono_spt(1,1)=Kjj;
      ext_mat_mono_spt(0,1)=K12;
      ext_mat_mono_spt(1,0)=K12;
      abs_vec_mono_spt[1] = interp(itw,scat_data.abs_vec_data(0,joker,joker,0,1),
                            t_gp,za_gp);

      if( stokes_dim == 2 ){
        break;
      }
      
      ext_mat_mono_spt(2,2)=Kjj;
      
      if( stokes_dim == 3 ){
        break;
      }
      
      K34=interp(itw,scat_data.ext_mat_data(0,joker,joker,0,2),t_gp,za_gp);
      ext_mat_mono_spt(2,3)=K34;
      ext_mat_mono_spt(3,2)=-K34;
      ext_mat_mono_spt(3,3)=Kjj;
      break;

    }
  default:
    out0 << "Not all particle type cases are implemented\n";
    
  }


}







//! pha_mat_singleCalc
/*!
 Returns the total phase matrix for given incident and scattered directions
. It requires a vector of particle number densities to be precalculated

 \param[out] Z               Output: phase matrix
 \param[out] za_sca          scattered 
 \param[out] aa_sca          and
 \param[out] za_inc          incident
 \param[out] aa_inc          directions
 \param[in]  scat_data_mono  workspace variable
 \param[in]  stokes_dim      workspace variable
 \param[in]  pnd_vec         vector of particle number densities at the point 
                             in question
 \param[in]  rte_temperature workspace variable
 \author Cory Davis
 \date   2003-11-27
*/
void pha_mat_singleCalc(
                        MatrixView& Z,                  
                        Numeric za_sca, 
                        Numeric aa_sca, 
                        Numeric za_inc, 
                        Numeric aa_inc,
                        const ArrayOfSingleScatteringData& scat_data_mono,
                        const Index&          stokes_dim,
                        const VectorView& pnd_vec,
                        const Numeric& rte_temperature
                        )
{
  Index N_pt=pnd_vec.nelem();

  assert(aa_inc>=-180 && aa_inc<=180);
  assert(aa_sca>=-180 && aa_sca<=180);

  Matrix Z_spt(stokes_dim, stokes_dim, 0.);
  Z=0.0;
  // this is a loop over the different particle types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      if (pnd_vec[i_pt]>0)
        {
          pha_mat_singleExtract(Z_spt,scat_data_mono[i_pt],za_sca,aa_sca,za_inc,
                                aa_inc,rte_temperature,stokes_dim);
          Z_spt*=pnd_vec[i_pt];
          Z+=Z_spt;
        }
    }
}

//! Extract the phase matrix from a monochromatic SingleScatteringData object
/*!
  Given a monochromatic SingleScatteringData object, incident and 
  scattered directions, and the temperature, this function returns the phase
  matrix in the laboratory frame
  \param[out] Z_spt the phase matrix
  \param[in]  scat_data a monochromatic SingleScatteringData object
  \param[in]  za_sca 
  \param[in]  aa_sca
  \param[in]  za_inc
  \param[in]  aa_inc
  \param[in]  rte_temperature
  \param[in]  stokes_dim

  \author Cory Davis
  \date 2004-07-16

*/
void pha_mat_singleExtract(
                           MatrixView& Z_spt,
                           const SingleScatteringData& scat_data,
                           const Numeric& za_sca,
                           const Numeric& aa_sca,
                           const Numeric& za_inc,
                           const Numeric& aa_inc,
                           const Numeric& rte_temperature,
                           const Index& stokes_dim
                           )                       
{
  switch (scat_data.ptype){

  case PTYPE_GENERAL:
    // to remove warnings during compilation. 
    out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
    break;
    
  case PTYPE_MACROS_ISO:
    {
      // Calculate the scattering and interpolate the data on the scattering
      // angle:
      
      Vector pha_mat_int(6);
      Numeric theta_rad;
            
      // Interpolation of the data on the scattering angle:
      interp_scat_angle_temperature(pha_mat_int, theta_rad, 
                                    scat_data, za_sca, aa_sca, 
                                    za_inc, aa_inc, rte_temperature);
      
      // Caclulate the phase matrix in the laboratory frame:
      pha_mat_labCalc(Z_spt, pha_mat_int, za_sca, aa_sca, za_inc, aa_inc, theta_rad);
      
      break;
    }

  case PTYPE_HORIZ_AL://Added by Cory Davis
    //Data is already stored in the laboratory frame, but it is compressed
    //a little.  Details elsewhere
    {
      assert (scat_data.pha_mat_data.ncols()==16);
      Numeric delta_aa=aa_sca-aa_inc+(aa_sca-aa_inc<-180)*360-
        (aa_sca-aa_inc>180)*360;//delta_aa corresponds to the "pages" 
                                //dimension of pha_mat_data
      GridPos t_gp;
      GridPos za_sca_gp;
      GridPos delta_aa_gp;
      GridPos za_inc_gp;
      Vector itw(16);
      gridpos(t_gp, scat_data.T_grid, rte_temperature);
      gridpos(delta_aa_gp,scat_data.aa_grid,abs(delta_aa));
      if (za_inc>90)
        {
          gridpos(za_inc_gp,scat_data.za_grid,180-za_inc);
          gridpos(za_sca_gp,scat_data.za_grid,180-za_sca);
        }
      else
        {
          gridpos(za_inc_gp,scat_data.za_grid,za_inc);
          gridpos(za_sca_gp,scat_data.za_grid,za_sca);
        }

      interpweights(itw,t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);

      Z_spt(0,0)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,0),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      if( stokes_dim == 1 ){
        break;
      }
      Z_spt(0,1)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,1),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      Z_spt(1,0)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,4),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      Z_spt(1,1)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,5),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      if( stokes_dim == 2 ){
        break;
      }
      if ((za_inc<=90 && delta_aa>=0)||(za_inc>90 && delta_aa<0))
        {
          Z_spt(0,2)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,2),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,2)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,6),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(2,0)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,8),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(2,1)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,9),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      else
        {
          Z_spt(0,2)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,2),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,2)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,6),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(2,0)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,8),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(2,1)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,9),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }                             
      Z_spt(2,2)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,10),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      if( stokes_dim == 3 ){
        break;
      }
      if ((za_inc<=90 && delta_aa>=0)||(za_inc>90 && delta_aa<0))
        {
          Z_spt(0,3)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,3),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,3)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,7),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,0)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,12),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,1)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,13),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      else
        {
          Z_spt(0,3)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,3),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(1,3)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,7),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,0)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,12),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
          Z_spt(3,1)=-interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                                   Range(joker),0,13),
                                  t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
        }
      Z_spt(2,3)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,11),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      Z_spt(3,2)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,14),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      Z_spt(3,3)=interp(itw,scat_data.pha_mat_data(0,joker,Range(joker),Range(joker),
                                               Range(joker),0,15),
                              t_gp,za_sca_gp,delta_aa_gp,za_inc_gp);
      break;
      
    }  
  default:
    out0 << "Not all particle type cases are implemented\n";
    
  }
}


//! Sample_los

/*!
  Implementation one of two line of sight sampling methods determined by the 
  input Index 'sampling_method'
  
  sampling_method==1:
     Randomly samples incident direction using a probability density function 
   proportional to sin(za)(cos(za)+1). sin(za) because dsolid_angle=sin(za)dzadaa
   , and (cos(za)+1) because upward radiances tend to be greater than downward
   radiances.  NOTE: THIS IS TERRIBLE IN OPTICALLY THICK CASES AND WILL BE 
   REMOVED

  sampling_method==2:
     Randomly samples incident direction using a probability density function 
   proportional to sin(za).  Better, but again, not very good.

  THIS WHOLE FUNCTION MAY BE SOON REDUNDANT TO SAMPLE_LOS_Z WHICH USES A PDF 
  PROPORTIONAL TO Z11SINZA

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
   \param[in]     rte_temperature

   \author Cory Davis
   \date   2003-06-19
*/

void Sample_los (
                   VectorView& new_rte_los,
                   Numeric& g_los_csc_theta,
                   MatrixView& Z,
                   Rng& rng,
                   const VectorView& rte_los,
                   const ArrayOfSingleScatteringData& scat_data_mono,
                   const Index&          stokes_dim,
                   const VectorView& pnd_vec,
                   const bool& anyptype30,
                   const VectorView& Z11maxvector,
                   Numeric Csca,
                   const Numeric& rte_temperature
                   )
{
  Numeric Z11max=0;
  bool tryagain=true;
  Numeric aa_scat = (rte_los[1]>=0) ?-180+rte_los[1]:180+rte_los[1];
      
  if(anyptype30)
    {
      Index np=pnd_vec.nelem();
      assert(scat_data_mono.nelem()==np);
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
      pha_mat_singleCalc(dummyZ,180-rte_los[0],aa_scat,180-rte_los[0],
                         aa_scat,scat_data_mono,stokes_dim,pnd_vec,rte_temperature);
      Z11max=dummyZ(0,0);
    }  
  ///////////////////////////////////////////////////////////////////////  
  while(tryagain)
    {
      new_rte_los[1] = rng.draw()*360-180;
      new_rte_los[0] = acos(1-2*rng.draw())*RAD2DEG;
      //Calculate Phase matrix////////////////////////////////
      Numeric aa_inc= (new_rte_los[1]>=0) ?
        -180+new_rte_los[1]:180+new_rte_los[1];
      
      pha_mat_singleCalc(Z,180-rte_los[0],aa_scat,180-new_rte_los[0],
                         aa_inc,scat_data_mono,stokes_dim,pnd_vec,rte_temperature);
      
      if (rng.draw()<=Z(0,0)/Z11max)//then new los is accepted
        {
          tryagain=false;
        }
    }
  g_los_csc_theta =Z(0,0)/Csca;
}


//! Sample_ppathlengthLOS

/*!
   Similar to Sample_ppathlength, but ensures that the sampled point lies within the cloudbox

  \param[out]    pathlength  the pathlength.
  \param[out]    g           the probability density of the returned pathlength.
  \param[in,out] rng         Rng random number generator instance
  \param[in]     TArray      FIXME: DOC
  \param[in]     cum_l_step  FIXME: DOC

  \author Cory Davis
  \date   2003-03-10
*/

void Sample_ppathlengthLOS (
                         Numeric& pathlength, 
                         Numeric& g,
                         Rng& rng,
                         const ArrayOfMatrix& TArray,
                         const ConstVectorView& cum_l_step
                         )
{
  Index npoints=cum_l_step.nelem();
  assert(TArray.nelem()==npoints);
        
  
  Numeric r = rng.draw();
  //inefficient first effort
  Vector T11vector(npoints);
  for (Index i=0;i<npoints;i++)
    {
      T11vector[i]=TArray[i](0,0);
      if (T11vector[i] < 1e-12)
        {
          //then we need to truncate T11vector at this point
          T11vector=T11vector[Range(0,i+1)];
          npoints=i+1;
          break;
        }
    }
  GridPos gp;
  //Vector itw(2);
  gridpos(gp,T11vector,1-r*(1-T11vector[npoints-1]));
  //interpweights(itw,gp);
  Numeric T11,T11A,T11B,sA,sB,K;
  T11=1-r*(1-T11vector[npoints-1]);
  T11A=T11vector[gp.idx];
  T11B=T11vector[gp.idx+1];
  sA=cum_l_step[gp.idx];
  sB=cum_l_step[gp.idx+1];
  K=log(T11A/T11B)/(sB-sA);//K=K11 only for diagonal ext_mat
  assert (K>0);
  pathlength=sA+log(T11A/T11)/K;
  g=K*T11/(1-T11vector[npoints-1]);
}               

//!  TArrayCalc
/*! 
 
   This function calculates arrays of useful variables along a propagation
   path within the cloudbox. 


   \author Cory Davis
   \date   2003-07-19
*/


                
void TArrayCalc(
                //output
                ArrayOfMatrix& TArray,
                ArrayOfMatrix& ext_matArray,
                ArrayOfVector& abs_vecArray,
                Vector& t_ppath,
                Tensor3& ext_mat,
                Matrix& abs_vec,
                Numeric&   rte_pressure,
                Numeric&   rte_temperature,
                Vector&    rte_vmr_list,
                Matrix&  pnd_ppath,
                //input
                const Ppath& ppath,
                const Agenda& opt_prop_gas_agenda,
                const Agenda& abs_scalar_gas_agenda,
                const Index& f_index,
                const Index& stokes_dim,
                const ConstVectorView&    p_grid_cloud,
                const ConstVectorView&    lat_grid_cloud,
                const ConstVectorView&    lon_grid_cloud,
                const ConstTensor3View&   t_field_cloud,
                const ConstTensor4View&   vmr_field_cloud,
                const Tensor4&   pnd_field,
                const ArrayOfSingleScatteringData& scat_data_mono,
                const ArrayOfIndex& cloudbox_limits
                )
{ 
  const Index np=ppath.np;  
  const Index   ns = vmr_field_cloud.nbooks();
  const Index N_pt = pnd_field.nbooks();
  
  Vector scat_za_grid(np);
  Vector scat_aa_grid(np);
  Vector p_ppath(np);
  Matrix vmr_ppath(ns,np);

  TArray.resize(np);
  ext_matArray.resize(np); 
  abs_vecArray.resize(np);
  pnd_ppath.resize(N_pt,np);
  t_ppath.resize(np);
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim,0.0);
  Matrix zeroMatrix(stokes_dim,stokes_dim,0.0);
  Matrix identity(stokes_dim,stokes_dim,0.0);
  //Identity matrix
  for (Index i=0; i<stokes_dim; i++){identity(i,i)=1.0;}
 
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);
  Vector abs_vec_part(stokes_dim, 0.0);

  //use propagation angles from ppath to form scat_za_grid, scat_aa_grid

  scat_za_grid=ppath.los(Range(joker),0);
  scat_za_grid*=-1.0;
  scat_za_grid+=180;
  scat_aa_grid=ppath.los(Range(joker),1);
  scat_aa_grid+=180;

  // Define cloud sub_grids and fields- this can actually happen one level above !!! ...
  
  cloud_atm_vars_by_gp(p_ppath,t_ppath,vmr_ppath,pnd_ppath,ppath.gp_p,
                       ppath.gp_lat,ppath.gp_lon,cloudbox_limits,p_grid_cloud,
                       lat_grid_cloud,lon_grid_cloud,t_field_cloud,
                       vmr_field_cloud,pnd_field);

  //Create array of extinction matrices corresponding to each point in ppath
  for (Index scat_za_index=0; scat_za_index<ppath.np; scat_za_index++)
    { 
      Matrix abs_scalar_gas;
      Index scat_aa_index=scat_za_index;
      rte_pressure    = p_ppath[scat_za_index];
      rte_temperature = t_ppath[scat_za_index];
      rte_vmr_list    = vmr_ppath(joker,scat_za_index);
      abs_scalar_gas_agendaExecute (abs_scalar_gas, f_index,
                                           rte_pressure, rte_temperature,
                                           rte_vmr_list,
                                           abs_scalar_gas_agenda, true);
      opt_prop_gas_agendaExecute( ext_mat, abs_vec, abs_scalar_gas,
                                  opt_prop_gas_agenda, true );
      ext_mat_part=0.0;
      abs_vec_part=0.0;
      //Make sure scat_aa is between -180 and 180
      if (scat_aa_grid[scat_aa_index]>180){scat_aa_grid[scat_aa_index]-=360;}
      //
      //opt_prop_part_agenda.execute( true );
      //use pnd_ppath and ext_mat_spt to get extmat (and similar for abs_vec
      
      opt_propCalc(ext_mat_part,abs_vec_part,scat_za_grid[scat_za_index],
                   scat_aa_grid[scat_aa_index],scat_data_mono,stokes_dim,
                   pnd_ppath(joker, scat_za_index),rte_temperature);

      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;
      abs_vec(0,Range(joker)) += abs_vec_part;

      ext_matArray[scat_za_index]=ext_mat(0,joker,joker);
      abs_vecArray[scat_za_index]=abs_vec(0,joker);
    }
  //create an array of T matrices corresponding to each position in the
  //the first ppath through the cloudbox each grid cell spanned by ppath points
  TArray[0]=identity;
  for (Index i=1; i<ppath.np;i++)
    {
      //Get the appropriate extinction matrix for the cell in question
      opt_depth_mat=ext_matArray[i];
      opt_depth_mat+=ext_matArray[i-1];
      opt_depth_mat*=-ppath.l_step[i-1]/2;
      incT=0;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index j=0;j<stokes_dim;j++)
            {
              incT(j,j)=exp(opt_depth_mat(j,j));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      TArray[i]=zeroMatrix;
      mult(TArray[i],TArray[i-1],incT);
    }

}


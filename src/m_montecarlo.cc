
/* Copyright (C) 2003 Cory Davis <cory@met.ed.ac.uk>
                            
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
  \file   m_montecarlo.cc
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2003-06-19 

  \brief  Workspace functions for the solution of cloud-box radiative transfer 
by Monte Carlo methods.  All of these functions refer to 3D calculations

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/
/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "messages.h"
#include "arts.h"
#include "ppath.h"
#include "matpackI.h"
#include "special_interp.h"
#include "check_input.h"
#include <stdexcept>
#include <cmath>
#include "rte.h"
#include "lin_alg.h"
#include "auto_md.h"
#include "logic.h"
#include "physics_funcs.h"
#include "xml_io.h"
#include "montecarlo.h"
#include "rng.h"
#include <ctime>
#include <fstream>
#include "mc_interp.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;
extern const Numeric BOLTZMAN_CONST;
extern const Numeric SPEED_OF_LIGHT;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! mc_errorApplySensor
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2005-03-15
*/
void mc_errorApplySensor(
           Vector&   mc_error,
     const Sparse&   sensor_response )
{
  const Index   n = mc_error.nelem();

  if( sensor_response.ncols() != n )
    {
      throw runtime_error(
                     "Mismatch in size of *sensor_response* and *mc_error*." );
    }

  Vector  i( n );
  for( Index j=0; j<n; j++ )
    { i[j] = mc_error[j]*mc_error[j]; }
  mc_error.resize( sensor_response.nrows() );
  mc_error = 0.0;
  for( Index irow=0; irow<sensor_response.nrows(); irow++ )
    {
      for( Index icol=0; icol<n; icol++ )
        { mc_error[irow] += pow( sensor_response(irow,icol), 2.0 ) * i[icol]; }
    }
  transform( mc_error, sqrt, mc_error );
}

void mc_antennaSetGaussian(
                           MCAntenna& mc_antenna,
                           //keyword arguments
                           const Numeric& za_sigma,
                           const Numeric& aa_sigma
                           )
{
  mc_antenna.set_gaussian(za_sigma,aa_sigma);
}

void mc_antennaSetGaussianByFWHM(
                           MCAntenna& mc_antenna,
                           //keyword arguments
                           const Numeric& za_fwhm,
                           const Numeric& aa_fwhm
                           )
{
  mc_antenna.set_gaussian_fwhm(za_fwhm,aa_fwhm);
}



void mc_antennaSetPencilBeam(
                            MCAntenna& mc_antenna
                            )
{
  mc_antenna.set_pencil_beam();
}



//! MCGeneral
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Cory Davis
   \date   2005-06-29
*/

void MCGeneral(
               Vector&               y,
               Index&                mc_iteration_count,
               Vector&               mc_error,
               Tensor3&              mc_points,
               MCAntenna&            mc_antenna,
               const Vector&         f_grid,
               const Matrix&         sensor_pos,
               const Matrix&         sensor_los,
               const Index&          stokes_dim,
               const Agenda&         iy_space_agenda,
               const Agenda&         surface_prop_agenda,
               const Agenda&         opt_prop_gas_agenda,
               const Agenda&         scalar_gas_absorption_agenda, 
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
               const Index&          mc_seed,
               const String&         mc_unit,
               //Keyword params
               const Numeric&        std_err,
               const Index&          max_time,
               const Index&          max_iter,
               const Index&          z_field_is_1D
               )
{
  //Check keyword input
  if (max_time<0 && max_iter<0 && std_err<0){
    throw runtime_error( "At least one of std_err, max_time, and max_iter must be positive" );
  }
  Ppath ppath_step;
  Rng rng;                      //Random Nuimber generator
  time_t start_time=time(NULL);
  Index N_pt=pnd_field.nbooks();//Number of particle types
  Vector pnd_vec(N_pt); //Vector of particle number densities used at each point
  Vector Z11maxvector;//Vector holding the maximum phase function for each 
  bool anyptype30=is_anyptype30(scat_data_mono);
  if (anyptype30)
    {
      findZ11max(Z11maxvector,scat_data_mono);
    }
  rng.seed(mc_seed);
  bool keepgoing,inside_cloud; // flag indicating whether to stop tracing a photons path
  Numeric g,temperature,albedo,g_los_csc_theta;
  Matrix A(stokes_dim,stokes_dim),Q(stokes_dim,stokes_dim),evol_op(stokes_dim,stokes_dim),
    ext_mat_mono(stokes_dim,stokes_dim),q(stokes_dim,stokes_dim),newQ(stokes_dim,stokes_dim),
    Z(stokes_dim,stokes_dim);
  q=0.0;newQ=0.0;
  mc_iteration_count=0;
  Vector vector1(stokes_dim),abs_vec_mono(stokes_dim),I_i(stokes_dim),Isum(stokes_dim),
    Isquaredsum(stokes_dim);
  y.resize(stokes_dim);
  y=0;
  Index termination_flag=0;
  mc_error.resize(stokes_dim);
  //local versions of workspace
  Matrix local_iy(1,stokes_dim),local_surface_emission(1,stokes_dim),local_surface_los;
  Tensor4 local_surface_rmatrix;
  Vector local_rte_pos(2);
  Vector local_rte_los(2);
  Vector new_rte_los(2);
    Index np;
  mc_points.resize(p_grid.nelem(),lat_grid.nelem(),lon_grid.nelem());
  mc_points=0;
  Isum=0.0;Isquaredsum=0.0;
  Numeric std_err_i;
  bool convert_to_rjbt=false;
  if ( mc_unit == "RJBT" )
    { 
      std_err_i=f_grid[0]*f_grid[0]*2*BOLTZMAN_CONST/SPEED_OF_LIGHT/SPEED_OF_LIGHT*
        std_err;
      convert_to_rjbt=true;
    }
  else if ( mc_unit == "radiance" )
    {
      std_err_i=std_err;
    }
  else
    {
      ostringstream os;
      os << "invalid value for mc_units - " << mc_unit <<". Must be \"RJBT\" or \"radiance\".";
      throw runtime_error( os.str() );
    }
      
  
  //Begin Main Loop
  while (true)
    {
      mc_iteration_count+=1;
      keepgoing=true; //flag indicating whether to continue tracing a photon path
      //Sample a FOV direction
      mc_antenna.draw_los(local_rte_los,rng,sensor_los(0,joker));
      id_mat(Q);
      local_rte_pos=sensor_pos(0,joker);
      I_i=0.0;
      while (keepgoing)
        {
          mcPathTraceGeneral(evol_op, abs_vec_mono, temperature, ext_mat_mono, rng, local_rte_pos, 
                             local_rte_los, pnd_vec, g,ppath_step,termination_flag, inside_cloud, 
                             opt_prop_gas_agenda,scalar_gas_absorption_agenda, 
                             stokes_dim, p_grid, 
                             lat_grid, lon_grid, z_field, r_geoid, z_surface,
                             t_field, vmr_field, cloudbox_limits, pnd_field,
                             scat_data_mono,z_field_is_1D);
          np=ppath_step.np;
          mc_points(ppath_step.gp_p[np-1].idx,ppath_step.gp_lat[np-1].idx,
                    ppath_step.gp_lon[np-1].idx)+=1;
          if (termination_flag==1)
            {
              iy_space_agendaExecute(local_iy,iy_space_agenda,true);
              mult(vector1,evol_op,local_iy(0,joker));
              mult(I_i,Q,vector1);
              I_i/=g;
              keepgoing=false; //stop here. New photon.
            }
          else if (termination_flag==2)
            {
              //decide whether we have reflection or emission
              surface_prop_agendaExecute(local_surface_emission, local_surface_los, 
                                         local_surface_rmatrix, ppath_step.gp_p[np-1],
                                         ppath_step.gp_lat[np-1],ppath_step.gp_lon[np-1],local_rte_los,
                                         surface_prop_agenda, true);
              //deal with blackbody case
              if (local_surface_los.nrows()==0)
                {
                  mult(vector1,evol_op,local_surface_emission(0,joker));
                  mult(I_i,Q,vector1);
                  I_i/=g;
                  keepgoing=false;
                }
              else
                //decide between reflection and emission
                {
                  Numeric R11=local_surface_rmatrix(0,0,0,0);
                  if (rng.draw()>R11)
                    {
                      //then we have emission
                      //Matrix oneminusR(stokes_dim,stokes_dim);
                      //id_mat(oneminusR);
                      //oneminusR-=local_surface_rmatrix(0,0,joker,joker);
                      //oneminusR/=1-R11;
                      //mult(vector1,oneminusR,local_surface_emission(0,joker));
                      mult(vector1,evol_op,local_surface_emission(0,joker));
                      mult(I_i,Q,vector1);
                      I_i/=g*(1-R11);
                      keepgoing=false;
                    }
                  else
                    {
                      //we have reflection
                      local_rte_los=local_surface_los(0,joker);
                      
                      mult(q,evol_op,local_surface_rmatrix(0,0,joker,joker));
                      mult(newQ,Q,q);
                      Q=newQ;
                      Q/=g*R11;
                    }
                }
            }
          else if (inside_cloud)
            {
              //we have another scattering/emission point 
              //Estimate single scattering albedo
              albedo=1-abs_vec_mono[0]/ext_mat_mono(0,0);
              //cout<<"albedo = "<<albedo<<" ext_mat_mono(0,0) = "<<ext_mat_mono(0,0)<<" abs_vec_mono[0] = "<<abs_vec_mono[0]<<"\n";
              //determine whether photon is emitted or scattered
              if (rng.draw()>albedo)
                {
                  //Calculate emission
                  Numeric planck_value = planck( f_grid[0], temperature );
                  Vector emission=abs_vec_mono;
                  emission*=planck_value;
                  Vector emissioncontri(stokes_dim);
                  mult(emissioncontri,evol_op,emission);
                  emissioncontri/=(g*(1-albedo));//yuck!
                  mult(I_i,Q,emissioncontri);
                  keepgoing=false;
                  //cout << "emission contri" <<  I_i[0] << "\n";
                }
              else
                {
                  //we have a scattering event
                  //Sample new line of sight.
                  
                  Sample_los (new_rte_los,g_los_csc_theta,Z,rng,local_rte_los,
                              scat_data_mono,stokes_dim,
                              pnd_vec,anyptype30,Z11maxvector,ext_mat_mono(0,0)-abs_vec_mono[0],temperature);
                                           
                  Z/=g*g_los_csc_theta*albedo;
                  
                  mult(q,evol_op,Z);
                  mult(newQ,Q,q);
                  Q=newQ;
                  //scattering_order+=1;
                  local_rte_los=new_rte_los;
                  //if (silent==0){cout <<"mc_iteration_count = "<<mc_iteration_count << 
                  //                 ", scattering_order = " <<scattering_order <<"\n";}
                }
            }
          else
            {
              //Must be clear sky emission point
              //Calculate emission
              Numeric planck_value = planck( f_grid[0], temperature );
              Vector emission=abs_vec_mono;
              emission*=planck_value;
              Vector emissioncontri(stokes_dim);
              mult(emissioncontri,evol_op,emission);
              emissioncontri/=g;
              mult(I_i,Q,emissioncontri);
              keepgoing=false;
              //cout << "emission contri" <<  I_i[0] << "\n";
            }
        }
      Isum += I_i;
      for(Index j=0; j<stokes_dim; j++)
        {
          assert(!isnan(I_i[j]));
          Isquaredsum[j] += I_i[j]*I_i[j];
        }
      y=Isum;
      y/=mc_iteration_count;
      for(Index j=0; j<stokes_dim; j++) 
        {
          mc_error[j]=sqrt((Isquaredsum[j]/mc_iteration_count-y[j]*y[j])/mc_iteration_count);
        }
      if (std_err>0 && mc_iteration_count>=100 && mc_error[0]<std_err_i){break;}
      if (max_time>0 && (Index)(time(NULL)-start_time)>=max_time){break;}
      if (max_iter>0 && mc_iteration_count>=max_iter){break;}
    }
  if ( convert_to_rjbt )
    {
      for(Index j=0; j<stokes_dim; j++) 
        {
          y[j]=invrayjean(y[j],f_grid[0]);
          mc_error[j]=invrayjean(mc_error[j],f_grid[0]);
        }
    }
}

//! MCSetSeedFromTime
/*!

    Sets the value of mc_seed from system time

   \author Cory Davis
   \date   2003-06-19
*/

void MCSetSeedFromTime(
                       Index&    mc_seed
                       )
{
  mc_seed=(Index)time(NULL);
}


//! ScatteringMonteCarlo
/*! 
 
   This workspace method uses a Monte Carlo method to calculate the monochromatic
   radiance leaving the cloudbox in one direction.  
   The theoretical basis for this algorithm can be found at 
   http://www.met.ed.ac.uk/~cory/davis_MRad04.pdf
   
   See the online help 
   (arts -d FUNCTION_NAME) for a description of parameters.
 

   \author Cory Davis
   \date   2003-06-19
*/

void ScatteringMonteCarlo (
                           // WS Output:
                           Ppath&                ppath,
                           Ppath&                ppath_step,
                           //Vector&               ppath_p,
                           //Vector&               ppath_t,
                           //Matrix&               ppath_vmr,
                           Vector&               mc_error,
                           Index&                mc_iteration_count,
                           Vector&               rte_pos,
                           Vector&               rte_los,
                           //Stuff needed by RteCalc
                           Matrix&               iy,
                           Numeric&              rte_pressure,
                           Numeric&              rte_temperature,
                           Vector&               rte_vmr_list,
                           //Other Stuff
                           Tensor3&              ext_mat,
                           Matrix&               abs_vec,
                           SLIData2& mc_incoming,
                           // WS Input:
                           const Agenda&         ppath_step_agenda,
                           const Index&          atmosphere_dim,
                           const Vector&         p_grid,
                           const Vector&         lat_grid,
                           const Vector&         lon_grid,
                           const Tensor3&        z_field,
                           const Matrix&         r_geoid,
                           const Matrix&         z_surface,
                           const ArrayOfIndex&   cloudbox_limits,
                           const Index&          stokes_dim,
                           //Stuff needed by RteCalc
                           const Agenda&         rte_agenda,
                           const Agenda&         iy_space_agenda,
                           const Agenda&         iy_surface_agenda,
                           const Tensor3&        t_field,
                           const Vector&         f_grid,
                           //Stuff needed by TArrayCalc
                           const Agenda& opt_prop_gas_agenda,
                           const Agenda& scalar_gas_absorption_agenda,
                           const Tensor4&   vmr_field,
                           //Other Stuff
                           const ArrayOfSingleScatteringData& scat_data_mono,
                           const Tensor4& pnd_field,
                           const Index& mc_seed,
                           const Index& f_index,
                           // Control Parameters:
                           const Numeric& std_err,
                           const Index& max_time,
                           const Index& max_iter,
                           const Index& incoming_lookup,
                           const Index& z_field_is_1D
                           )

{ 
  /*These used to be keyword parameters, but they are so seldom used that they 
    are now hard coded. The next step would be to remove all of the code that 
    they refer to.*/
  const Index& record_ppathcloud=0;
  const Index& record_ppath=0;
  const Index& silent=1;
  const Index& record_histdata=0;
  const String& histdata_filename="";
                           
  //Check keyword input
  if (max_time<0 && max_iter<0 && std_err<0){
    throw runtime_error( "At least one of std_err, max_time, and max_iter must be positive" );
  }
              
  //INTERNAL DECLARATIONS/////////////////////////////////////////////////
  Matrix identity(stokes_dim,stokes_dim,0.0); //The identity matrix
  for (Index i=0; i<stokes_dim; i++){identity(i,i)=1.0;}

  Matrix Q(stokes_dim,stokes_dim);//defined in eq 12 of ref1
  Matrix evol_op(stokes_dim,stokes_dim);//evolution operator (symbol O in ref1)
  Matrix K(stokes_dim,stokes_dim);//extinction matrix
  bool keepgoing; // flag indicating whether to stop tracing a photons path
  Index scattering_order;       //k in ref1
  Vector new_rte_los(2);        //Stores new line of sight
  Ppath ppathLOS;               // propagation path in the original line of sight
  Ppath ppathcloud;             //prop. path inside cloud box
  Numeric g;             //Probability density for pathlength sampling
  Numeric pathlength;           //\Delta s in ref1
  ArrayOfMatrix TArrayLOS;//array of evolution operators along 
                                //original line of sight                
  ArrayOfMatrix TArray;//array of evolution operators for higher 
                                //scattering orders
  ArrayOfMatrix ext_matArray;//array of extinction matrices along 
                                //propagation path
  ArrayOfMatrix ext_matArrayLOS;//array of extinction matrices along
                                //original LOS
  ArrayOfVector abs_vecArray;//array of abs. coeff. vectors along 
                                //propagation path
  ArrayOfVector abs_vecArrayLOS;//array of abs. coeff. vectors along 
                                //original LOS
  Vector Isum(stokes_dim,0.0);//Sum of all photon contributions to the Stokes vector
  Vector Isquaredsum(stokes_dim,0.0);//Used to estimate error
  Vector Iboundary(stokes_dim,0.0);//Incoming Stokes vector at the cloudbox boundary
  Vector IboundaryLOScontri(stokes_dim,0.0);//1st term RHS eq 2 ref1.
  Matrix Z(stokes_dim,stokes_dim,0.0);//bulk phase matrix
  Matrix q(stokes_dim,stokes_dim,0.0);//Eqs. 12-13 ref1
  Matrix newQ(stokes_dim,stokes_dim,0.0);//Eq 12 ref1
  Vector cum_l_step;//vector of cumulative distance along a propagation path
  Vector cum_l_stepLOS;// "    "  for the original LOS
  Vector t_ppath;//vector of temperatures along a propagation path
  Vector t_ppathLOS;//" "  " along the original LOS
  Matrix pnd_ppath;//ppath.np by Nptypes matrix of pasrticle number density
  Matrix pnd_ppathLOS;// "   "    "
  ArrayOfGridPos pathlength_gp(1);//for optical properties along a propagation 
                                //path according to pathlength
  Vector K_abs(stokes_dim);//absorption coefficient vector
  Vector I(stokes_dim);//Cloudbox exit Stokes Vector
  mc_error.resize(stokes_dim);//Error in Cloudbox exit Stokes vector
  Rng rng;                      //Random Nuimber generator
  Vector I_i(stokes_dim);//photon contribution to Stokes vector
  Vector boundarycontri(stokes_dim);//eq 16 ref1
  Numeric g_los_csc_theta;//eq. 11 ref1 divided by sin(\theta) 
  Numeric albedo;//eq. 9 ref1
  Numeric dist_to_boundary; //Distance to the far boundary of the cloudbox
  Index N_pt=pnd_field.nbooks();//Number of particle types
  Vector pnd_vec(N_pt); //Vector of particle number densities used at each point
  time_t start_time=time(NULL);
  Vector Z11maxvector;//Vector holding the maximum phase function for each 
  //particle type. Used in Rejection method for sampling incident direction
  bool left_cloudbox;
  //////////////////////////////////////////////////////////////////////////////
  //Calculate the clea-sky transmittance between cloud box and sensor using the
  //existing ppath.
  Numeric transmittance=
    exp(-opt_depth_calc(ext_mat, abs_vec, rte_pressure, rte_temperature, 
                        rte_vmr_list, ppath, opt_prop_gas_agenda,
                        scalar_gas_absorption_agenda, f_index, p_grid,
                        lat_grid, lon_grid, t_field, vmr_field,
                        atmosphere_dim));
  Numeric std_err_i=f_grid[0]*f_grid[0]*2*BOLTZMAN_CONST/SPEED_OF_LIGHT/SPEED_OF_LIGHT*
    std_err/transmittance;
                            
  bool anyptype30=is_anyptype30(scat_data_mono);
  if (anyptype30)
    {
      findZ11max(Z11maxvector,scat_data_mono);
    }
  //If necessary, open file for histogram data output
  ofstream histfile;
  if (record_histdata==1)
    {
      const char* p = histdata_filename.c_str();
      histfile.open(p,ios::out);
    }
  //out1 << rte_pos << "\n";
  //out1 << rte_los << "\n";
  rng.seed(mc_seed);
  Agenda iy_cloudbox_agenda;
  Cloudbox_ppath_rteCalc(ppathLOS, ppath, ppath_step, 
                         //ppath_p, ppath_t, ppath_vmr, 
                         rte_pos, rte_los, 
                         cum_l_stepLOS, TArrayLOS, ext_matArrayLOS, 
                         abs_vecArrayLOS,t_ppathLOS, ext_mat, abs_vec, rte_pressure, 
                         rte_temperature, rte_vmr_list, iy, pnd_ppathLOS, 
                         ppath_step_agenda, atmosphere_dim, 
                         p_grid, lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                         cloudbox_limits, record_ppathcloud, record_ppath, 
                         opt_prop_gas_agenda, 
                         scalar_gas_absorption_agenda, f_index,
                         stokes_dim, t_field, vmr_field, rte_agenda,
                         iy_space_agenda, 
                         iy_surface_agenda, iy_cloudbox_agenda,f_grid, 0, 0,
                         pnd_field,scat_data_mono, z_field_is_1D );
 
  mult(IboundaryLOScontri,TArrayLOS[TArrayLOS.nelem()-1],iy(0,joker));

  //out1 << iy << "\n";
  //out1 << TArrayLOS[TArrayLOS.nelem()-1] << "\n";
  //out1 << TArrayLOS.nelem() << "\n";

  for (Index i = 0;i<stokes_dim;i++){assert(!isnan(IboundaryLOScontri[i]));}
  
  mc_iteration_count=0;
  //Begin Main Loop
  while (true)
    {
      mc_iteration_count+=1;
      left_cloudbox=false;
      keepgoing=true; //flag indicating whether to continue tracing a photon path
      scattering_order=0;       //scattering order
      Q=identity;               //identity matrix
      I_i=0.0;
      boundarycontri=0.0;
      //while the reversed traced photon path remains in the cloud box
      //
      TArray=TArrayLOS;
      ext_matArray=ext_matArrayLOS;
      abs_vecArray=abs_vecArrayLOS;
      ppathcloud=ppathLOS;
      cum_l_step=cum_l_stepLOS;
      t_ppath=t_ppathLOS;
      pnd_ppath=pnd_ppathLOS;
      dist_to_boundary=cum_l_step[ppathcloud.np-1];
         
      if (silent==0){cout<<"mc_iteration_count = "<<mc_iteration_count<<"\n";}
      while (keepgoing)
        {
           if (scattering_order>0)
            {
              mcPathTrace(evol_op, K_abs, rte_temperature, K, rng, rte_pos, 
                          rte_los, pnd_vec, g,left_cloudbox, opt_prop_gas_agenda,
                          scalar_gas_absorption_agenda, stokes_dim, p_grid, 
                          lat_grid, lon_grid, z_field, r_geoid, z_surface,
                          t_field, vmr_field, cloudbox_limits, pnd_field,
                          scat_data_mono,z_field_is_1D);
            }
          else
            {
              Sample_ppathlengthLOS (pathlength,g,rng,TArray,cum_l_step);
              //Interpolate T, s_0, etc from ppath and Tarray
              interpTArray(evol_op, K_abs, rte_temperature, K, rte_pos, rte_los, 
                           pnd_vec, TArray, ext_matArray,
                           abs_vecArray, t_ppath, pnd_ppath, cum_l_step,
                           pathlength, stokes_dim, ppathcloud);
            }
          assert(cum_l_step.nelem()==ppathcloud.np);
          assert(TArray.nelem()==ppathcloud.np);
          if ( left_cloudbox )
            //Then the path has left the cloud box
            {
              //Get incoming//////
              if ( incoming_lookup )
                {
                  Iboundary=0;
                  //Iboundary[0]=mcGetIncomingFromData(rte_pos,rte_los,mc_incoming);
                  Iboundary[0]=mc_incoming.interpolate(rte_pos[0],rte_los[0]);
                }
              else
                {
                  montecarloGetIncoming(iy,rte_pos,rte_los,
                                        ppath,ppath_step,
                                        //ppath_p, ppath_t, ppath_vmr,
                                        ppath_step_agenda,
                                        rte_agenda,iy_space_agenda,iy_surface_agenda,
                                        iy_cloudbox_agenda,
                                        p_grid,lat_grid,lon_grid,z_field,
                                        t_field, vmr_field, r_geoid,
                                        z_surface,cloudbox_limits, atmosphere_dim,
                                        f_grid,stokes_dim);
                  
                  Iboundary=iy(0,joker);
                }
              ////////////////////
              mult(boundarycontri,evol_op,Iboundary);
              mult(I_i,Q,boundarycontri);
              I_i/=g;
              //cout << "boundary contri" <<  I_i[0] << "\n";
              keepgoing=false; //stop here. New photon.
            }
          else
            {
              //we have another scattering/emission point 
              //Estimate single scattering albedo
              albedo=1-K_abs[0]/K(0,0);
              //cout<<"albedo = "<<albedo<<" K(0,0) = "<<K(0,0)<<" K_abs[0] = "<<K_abs[0]<<"\n";
              //determine whether photon is emitted or scattered
              if (rng.draw()>albedo)
                {
                  //Calculate emission
                  Numeric planck_value = planck( f_grid[f_index], rte_temperature );
                  Vector emission=K_abs;
                  emission*=planck_value;
                  Vector emissioncontri(stokes_dim);
                  mult(emissioncontri,evol_op,emission);
                  emissioncontri/=(g*(1-albedo));//yuck!
                  mult(I_i,Q,emissioncontri);
                  keepgoing=false;
                  //cout << "emission contri" <<  I_i[0] << "\n";
                }
              else
                {
                  //Sample new line of sight.
                  
                  Sample_los (new_rte_los,g_los_csc_theta,Z,rng,rte_los,
                              scat_data_mono,stokes_dim,
                              pnd_vec,anyptype30,Z11maxvector,K(0,0)-K_abs[0],rte_temperature);
                                           
                  Z/=g*g_los_csc_theta*albedo;
                  
                  mult(q,evol_op,Z);
                  mult(newQ,Q,q);
                  Q=newQ;
                  scattering_order+=1;
                  rte_los=new_rte_los;
                  if (silent==0){cout <<"mc_iteration_count = "<<mc_iteration_count << 
                                   ", scattering_order = " <<scattering_order <<"\n";}
                }

            }
 
        }
      Isum += I_i;
      if (record_histdata==1){histfile << I_i << "\n";}
      for(Index j=0; j<stokes_dim; j++)
        {
          assert(!isnan(I_i[j]));
          Isquaredsum[j] += I_i[j]*I_i[j];
        }
      I=Isum;
      I/=mc_iteration_count;
      for(Index j=0; j<stokes_dim; j++) 
        {
          mc_error[j]=sqrt((Isquaredsum[j]/mc_iteration_count-I[j]*I[j])/mc_iteration_count);
        }
      if (std_err>0 && mc_iteration_count>=100 && mc_error[0]<std_err_i){break;}
      if (max_time>0 && (Index)(time(NULL)-start_time)>=max_time){break;}
      if (max_iter>0 && mc_iteration_count>=max_iter){break;}
    }
  
  I+=IboundaryLOScontri;
  iy(0,joker)=I;
  mc_error*=transmittance;
}               




//! rte_posShift
/*! 
   shifts rte_pos and rte_los, and rte_gp_XXX to the end of ppath.

   \author Cory Davis
   \date   2003-07-19
*/


void rte_posShift(
                  Vector&         rte_pos,
                  Vector&         rte_los,
                  GridPos&        rte_gp_p,
                  GridPos&        rte_gp_lat,
                  GridPos&        rte_gp_lon,
                  const Ppath&    ppath,
                  const Index&    atmosphere_dim)
{
  const Index np      = ppath.np;
  
  rte_pos.resize( atmosphere_dim );
  rte_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  rte_los.resize( ppath.los.ncols() );
  rte_los = ppath.los(np-1,joker);
  gridpos_copy( rte_gp_p, ppath.gp_p[np-1] );
  if( atmosphere_dim > 1 )
    { gridpos_copy( rte_gp_lat, ppath.gp_lat[np-1] ); }
  if( atmosphere_dim > 2 )
    { gridpos_copy( rte_gp_lon, ppath.gp_lon[np-1] ); }
}  




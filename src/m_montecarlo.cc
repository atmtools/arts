
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

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;
extern const Numeric BOLTZMAN_CONST;
extern const Numeric SPEED_OF_LIGHT;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! scat_iPutMonteCarlo
/*! 
calculates interface Tensors scat_i_p, scat_i_lat, and scat_i_lon.  This is
the equivalent of scat_iPut for use after ScatteringMonteCarlo and before a final
call to RteCalc.  See the online help (arts -d FUNCTION_NAME) for a description 
of parameters.
   
\author Cory Davis
\date 2003-06-30
  
*/

void scat_iPutMonteCarlo(
                         Tensor7& scat_i_p,
                         Tensor7& scat_i_lat,
                         Tensor7& scat_i_lon,
                         const Matrix& iy,
                         const Index& stokes_dim,
                         const Vector& f_grid,
                         const ArrayOfIndex& cloudbox_limits,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid
                         )
{

  Vector I = iy(0,joker);
   Index Nf = f_grid.nelem();
   Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

   Index Nza = scat_za_grid.nelem();
   
   Index Ni = stokes_dim;
   Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1; 
   Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
   Index Naa = scat_aa_grid.nelem();
   
   scat_i_p.resize(Nf, 2, Nlat_cloud, Nlon_cloud, Nza, Naa, Ni);
   scat_i_lat.resize(Nf, Np_cloud, 2, Nlon_cloud, Nza, Naa, Ni);
   scat_i_lon.resize(Nf, Np_cloud, Nlat_cloud, 2, Nza, Naa, Ni);


  for(Index i = 0;i<stokes_dim;i++)
    {
      scat_i_p(joker,joker,joker,joker,joker,joker,i)=I[i];
      scat_i_lat(joker,joker,joker,joker,joker,joker,i)=I[i];
      scat_i_lon(joker,joker,joker,joker,joker,joker,i)=I[i];
    }
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
                           Vector&               mc_error,
                           Index&                mc_iteration_count,
                           Vector&               rte_pos,
                           Vector&               rte_los,
                           //Stuff needed by RteCalc
                           GridPos&              rte_gp_p,
                           GridPos&              rte_gp_lat,
                           GridPos&              rte_gp_lon,
                           // Matrix&               i_space,
                           //Matrix&               surface_emission,
                           //Matrix&               surface_los, 
                           //Tensor4&              surface_refl_coeffs,
                           Matrix&               iy,
                           //Vector&               scat_za_grid,
                           //Vector&               scat_aa_grid,
                           Numeric&              rte_pressure,
                           Numeric&              rte_temperature,
                           Vector&               rte_vmr_list,
                           //Other Stuff
                           Tensor3&              ext_mat,
                           Matrix&               abs_vec,
                           Index&                f_index,
                          
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
                           // Control Parameters:
                           const Numeric& std_err,
                           const Index& max_time,
                           const Index& max_iter,
                           const Index& rng_seed
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
  Matrix T(stokes_dim,stokes_dim);//evolution operator (symbol O in ref1)
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
  Numeric za_scat;//zenith angle of scattering direction
  Vector Z11maxvector;//Vector holding the maximum phase function for each 
  //particle type. Used in Rejection method for sampling incident direction
  //////////////////////////////////////////////////////////////////////////////
  //Calculate the clea-sky transmittance between cloud box and sensor using the
  //existing ppath.
  Numeric transmittance=exp(-opt_depth_calc(ext_mat, rte_pressure, rte_temperature, 
                                       rte_vmr_list, ppath, opt_prop_gas_agenda,
                                       scalar_gas_absorption_agenda, p_grid,
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
  //if rng_seed is < 0, keep time based seed, otherwise...
  if(rng_seed>=0){rng.seed(rng_seed);}
  Agenda iy_cloudbox_agenda;
  Cloudbox_ppath_rteCalc(ppathLOS, ppath, ppath_step, rte_pos, rte_los, 
                         cum_l_stepLOS, TArrayLOS, ext_matArrayLOS, 
                         abs_vecArrayLOS,t_ppathLOS, ext_mat, abs_vec, rte_pressure, 
                         rte_temperature, rte_vmr_list, iy, rte_gp_p, 
                         rte_gp_lat, rte_gp_lon, f_index, pnd_ppathLOS, 
                         ppath_step_agenda, atmosphere_dim, 
                         p_grid, lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                         cloudbox_limits, record_ppathcloud, record_ppath, 
                         opt_prop_gas_agenda, 
                         scalar_gas_absorption_agenda, stokes_dim, t_field, 
                         vmr_field, rte_agenda, iy_space_agenda, 
                         iy_surface_agenda, iy_cloudbox_agenda,f_grid, 0, 0,
                         pnd_field,scat_data_mono);
 
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
              //We need to calculate a new propagation path. In the future, we may be 
              //able to take some shortcuts here
              //Cloudbox_ppathCalc(ppathcloud,ppath_step,ppath_step_agenda,atmosphere_dim,
              //               p_grid,lat_grid,lon_grid,z_field,r_geoid,z_surface,
              //              cloudbox_limits, rte_pos,rte_los);
              ppath_calc(ppathcloud,ppath_step,ppath_step_agenda,
                      atmosphere_dim,p_grid,lat_grid,lon_grid,z_field,
                       r_geoid,z_surface,1,cloudbox_limits, rte_pos,rte_los,0);
          
              if (record_ppathcloud){ppathRecordMC(ppathcloud,"ppathcloud",
                                                   mc_iteration_count,scattering_order);}
                              
              cum_l_stepCalc(cum_l_step,ppathcloud);
  
              //Calculate array of transmittance matrices
              TArrayCalc(TArray, ext_matArray, abs_vecArray, t_ppath, ext_mat, abs_vec, 
                         rte_pressure, rte_temperature, 
                         rte_vmr_list, pnd_ppath, ppathcloud, opt_prop_gas_agenda, 
                         scalar_gas_absorption_agenda, stokes_dim, 
                         p_grid, lat_grid, lon_grid, t_field, vmr_field, atmosphere_dim,
                         pnd_field, scat_data_mono);


              /////////////////////////////////////////////////////////////////////
              dist_to_boundary=cum_l_step[ppathcloud.np-1];
         
              //              Iboundary=iy(0,joker);
              Sample_ppathlength (pathlength,g,rng,TArray,cum_l_step);
            }
          else
            {
              Sample_ppathlengthLOS (pathlength,g,rng,TArray,cum_l_step);
            }
          assert(cum_l_step.nelem()==ppathcloud.np);
          assert(TArray.nelem()==ppathcloud.np);
          if (pathlength>dist_to_boundary)
            //Then the path has left the cloud box
            {
              assert (scattering_order>0); //scattering/emission should be 
                                           //forced in original line of sight
              //Get incoming//////
              montecarloGetIncoming(iy,rte_pos,rte_los,rte_gp_p,
                        rte_gp_lat,rte_gp_lon,ppath,ppath_step,
                                    ppath_step_agenda,
                        rte_agenda,iy_space_agenda,iy_surface_agenda,
                                    iy_cloudbox_agenda,
                        p_grid,lat_grid,lon_grid,z_field,r_geoid,
                        z_surface,cloudbox_limits,ppathcloud,atmosphere_dim,
                        f_grid,stokes_dim);
              
              if (record_ppath)
                {
                  ppathRecordMC(ppath,"ppath",mc_iteration_count,scattering_order);
                }
              
              f_index=0;//For some strange reason f_index is set to -1 in RteStandard
              Iboundary=iy(0,joker);
              ////////////////////
              T=TArray[ppathcloud.np-1];
              mult(boundarycontri,T,Iboundary);
              mult(I_i,Q,boundarycontri);
              I_i/=g;
              keepgoing=false; //stop here. New photon.
            }
          else
            {
               //we have another scattering/emission point
              //Interpolate T, s_0, etc from ppath and Tarray
              interpTArray(T, K_abs, rte_temperature, K, rte_pos, rte_los, pnd_vec,
                           pathlength_gp,TArray, 
                           ext_matArray,abs_vecArray, t_ppath, pnd_ppath,
                           cum_l_step,pathlength, 
                           stokes_dim, ppathcloud);
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
                  mult(emissioncontri,T,emission);
                  emissioncontri/=(g*(1-albedo));//yuck!
                  mult(I_i,Q,emissioncontri);
                  keepgoing=false;
                   
                }
              else
                {
                  //Sample new line of sight.
                  
                  Sample_los (new_rte_los,g_los_csc_theta,Z,rng,rte_los,
                              scat_data_mono,stokes_dim,
                              pnd_vec,anyptype30,Z11maxvector,K(0,0)-K_abs[0],rte_temperature);
                                           
                  Z/=g*g_los_csc_theta*albedo;
                  
                  mult(q,T,Z);
                  mult(newQ,Q,q);
                  Q=newQ;
                  scattering_order+=1;
                  za_scat=180-rte_los[0];
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





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

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


//! Cloudbox_ppath_calc
/*! 
  This function performs the same task as ppath_calc, except inside the
  cloudbox.  It has been derived from the clear sky version.  See the 
  online help (arts -d FUNCTION_NAME) for a description of parameters.
   
\author Cory Davis
\date 2003-06-19
  
*/
void Cloudbox_ppathCalc(
        // WS Output:
              Ppath&          ppath,
              Ppath&          ppath_step,
        // WS Input:
        const Agenda&         ppath_step_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         rte_pos,
        const Vector&         rte_los)
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
      chk_if_in_range( "sensor latitude", rte_pos[1], -90, 90 );
      chk_if_in_range( "sensor longitude", rte_pos[2], -360, 360 );
      chk_vector_length( "rte_los", rte_los, 2 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], 0, 180 );
      chk_if_in_range( "sensor azimuth angle", rte_los[1], -180, 180 );
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
                      lon_grid, z_field, r_geoid, z_ground, rte_pos, rte_los );

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
      ppath_step_agenda.execute( true );

      // Before everything is tested carefully, we consider more than 1000
      // path points to be an indication on that the calcululations have
      // got stuck in an infinite loop.
      if( istep > 5000 )
        {
          throw logic_error(
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


//! montecarlo_p_from_belowCscaAdapt
/*! 
Reduces the montecarlo_p_from_belowCsca Tensor to the frequency of the 
Monte Carlo simulation.
   
\author Cory Davis
\date 2003-12-4
  
*/

void montecarlo_p_from_belowCscaAdapt(
				      Tensor3& montecarlo_p_from_belowCsca,
				      const Vector& f_grid,
				      const Index& f_index,
				      const ArrayOfSingleScatteringData& scat_data_raw
				      )
{
  Index N_pt=montecarlo_p_from_belowCsca.npages();
  //  N_f = montecarlo_p_from_belowCsca.ncols();
  Index N_za=montecarlo_p_from_belowCsca.ncols();
  Tensor3 new_p_from_belowCsca=Tensor3(N_pt,1,N_za);
  Numeric freq=f_grid[f_index];
  GridPos gp;
  Vector itw(2);
  Vector old_grid=scat_data_raw[0].f_grid;
  for(Index pt_index=0;pt_index<N_pt;pt_index++)
    {
      for (Index za_index=0;za_index<N_za;za_index++)
	{
	  gridpos(gp,old_grid,freq);
	  interpweights(itw,gp);
	  new_p_from_belowCsca(pt_index,0,za_index)=interp(itw,
		 montecarlo_p_from_belowCsca(pt_index,Range(joker),za_index),gp);
	}
    }
  montecarlo_p_from_belowCsca=new_p_from_belowCsca;
}









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
			 const Matrix& i_rte,
			 const Index& stokes_dim,
                         const Vector& f_grid,
                         const ArrayOfIndex& cloudbox_limits,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid
			 )
{

  Vector I = i_rte(0,joker);
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
   radiance leaving the cloudbox in one direction.  At this stage it is a very basic
   function - the phase matrix is not used to sample incident propagation directions.
   For phase matrix LOS sampling to be worthwhile, other features need to be added, 
   for example stratified sampling.  In the meantime, this function should give the 
   correct result, but may not be very efficient.See the online help 
   (arts -d FUNCTION_NAME) for a description of parameters.
 

   \author Cory Davis
   \date   2003-06-19
*/

void ScatteringMonteCarlo (
			   // WS Output:
			   Ppath&                ppath,
			   Ppath&                ppath_step,
			   Vector&               i_montecarlo_error,
			   Vector&               rte_pos,
			   Vector&               rte_los,
			   //Stuff needed by RteCalc
			   GridPos&              rte_gp_p,
			   GridPos&              rte_gp_lat,
			   GridPos&              rte_gp_lon,
			   Matrix&               i_space,
			   Matrix&               ground_emission,
			   Matrix&               ground_los, 
			   Tensor4&              ground_refl_coeffs,
			   Matrix&               i_rte,
			   Vector&               scat_za_grid,
			   Vector&               scat_aa_grid,
			   Numeric&              rte_pressure,
			   Numeric&              rte_temperature,
			   Vector&               rte_vmr_list,
			   //Other Stuff
			   Tensor3&              ext_mat,
			   Matrix&               abs_vec,
			   Index&                f_index,
			   Index&                scat_za_index,
			   Index&                scat_aa_index,
			   Tensor3&              ext_mat_spt,
			   Matrix&               abs_vec_spt,
			  
			   // WS Input:
			   const Agenda&         ppath_step_agenda,
			   const Index&          atmosphere_dim,
			   const Vector&         p_grid,
			   const Vector&         lat_grid,
			   const Vector&         lon_grid,
			   const Tensor3&        z_field,
			   const Matrix&         r_geoid,
			   const Matrix&         z_ground,
			   const ArrayOfIndex&   cloudbox_limits,
			   const Index&          stokes_dim,
			   //Stuff needed by RteCalc
			   const Agenda&         rte_agenda,
			   const Agenda&         i_space_agenda,
			   const Agenda&         ground_refl_agenda,
			   const Tensor3&        t_field,
			   const Vector&         f_grid,
			   //Stuff needed by TArrayCalc
			   const Agenda& opt_prop_gas_agenda,
			   const Agenda& spt_calc_agenda,
			   const Agenda& scalar_gas_absorption_agenda,
			   const Tensor4&   vmr_field,
                           //Other Stuff
			   const ArrayOfSingleScatteringData& scat_data_raw,
			   const Tensor4& pnd_field,
                           const Tensor4& scat_theta, // CE: Included 
                           const ArrayOfArrayOfArrayOfArrayOfGridPos& scat_theta_gps,
                           const Tensor5& scat_theta_itws,
			   const Tensor3& montecarlo_p_from_belowCsca,
			    // Control Parameters:
			   const Index& maxiter,
			   const Index& rng_seed,
			   const Index& record_ppathcloud,
			   const Index& record_ppath,
			   const Index& silent,
			   const Index& record_histdata,
			   const String& histdata_filename,
			   const Index& los_sampling_method,
			   const Index& strat_sampling
                           )

{		
  //Internal Declarations
   Matrix identity(stokes_dim,stokes_dim,0.0);
  //Identity matrix
  for (Index i=0; i<stokes_dim; i++){identity(i,i)=1.0;}
  Matrix Q(stokes_dim,stokes_dim),T(stokes_dim,stokes_dim);
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim);
  Matrix K(stokes_dim,stokes_dim);
  bool keepgoing;
  Index scattering_order;
  Index photon_number;
  Vector new_rte_los(2);
  Ppath ppathLOS,ppathcloud;
  Numeric g;
  Numeric pathlength;
  Numeric temperature;		
  ArrayOfMatrix TArrayLOS(ppath.np);
  ArrayOfMatrix TArray(ppath.np);
  ArrayOfMatrix ext_matArray(ppath.np);
  ArrayOfMatrix ext_matArrayLOS(ppath.np);
  ArrayOfVector abs_vecArray(ppath.np);
  ArrayOfVector abs_vecArrayLOS(ppath.np);
  Vector Isum(stokes_dim,0.0);
  Vector Isquaredsum(stokes_dim,0.0);
  Vector Iboundary(stokes_dim,0.0);
  Vector IboundaryLOScontri(stokes_dim,0.0);
  Vector pha_mat_za_grid(2);
  Vector pha_mat_aa_grid(2);
  Matrix Z(stokes_dim,stokes_dim,0.0);
  Matrix q(stokes_dim,stokes_dim,0.0);
  Matrix newQ(stokes_dim,stokes_dim,0.0);
  Vector cum_l_step;
  Vector cum_l_stepLOS;
  Vector t_ppath;
  Vector t_ppathLOS;
  Matrix pnd_ppath;
  Matrix pnd_ppathLOS;
  Tensor5 pha_mat_spt(scat_data_raw.nelem(),2,2,stokes_dim,stokes_dim);
  Tensor4 pha_mat(2,2,stokes_dim,stokes_dim);
  ArrayOfGridPos pathlength_gp(1);
  Vector K_abs(stokes_dim);
  Vector I(stokes_dim);
  i_montecarlo_error.resize(stokes_dim);
  Rng rng;
  Vector pathI(stokes_dim);
  Vector boundarycontri(stokes_dim);
  Vector pathinc(stokes_dim);
  Numeric g_los_csc_theta;
  Numeric albedo;
  Numeric dist_to_boundary;
  Numeric K11;
  Index N_pt=pnd_field.nbooks();
  Vector pnd_vec(N_pt);
  time_t start_time=time(NULL);
  //Stratified sampling variables//////////
  Numeric p_from_below;
  Vector I_from_below_sum(stokes_dim,0.0);
  Vector I_from_below_squaredsum(stokes_dim,0.0);
  Index I_from_below_N=0;
  Vector I_from_above_sum(stokes_dim,0.0);
  Vector I_from_above_squaredsum(stokes_dim,0.0);
  Index I_from_above_N=0;
  Vector I_emission_sum(stokes_dim,0.0);
  Vector I_emission_squaredsum(stokes_dim,0.0);
  Index I_emission_N=0;
  Numeric za_scat;
  Vector za_grid=scat_data_raw[0].za_grid;
  /////////////////////////////////////////

  //If necessary, open file for histogram data output
  ofstream histfile;
  if (record_histdata==1)
    {
      const char* p = histdata_filename.c_str();
      histfile.open(p,ios::out);
    }

  //if rng_seed is < 0, keep time based seed, otherwise...
  if(rng_seed>=0){rng.seed(rng_seed);}
  
  Cloudbox_ppath_rteCalc(ppathLOS, ppath, ppath_step, rte_pos, rte_los, 
			 cum_l_stepLOS, TArrayLOS, ext_matArrayLOS, 
			 abs_vecArrayLOS,t_ppathLOS, scat_za_grid, 
			 scat_aa_grid, ext_mat, abs_vec, rte_pressure, 
			 rte_temperature, rte_vmr_list, i_rte, rte_gp_p, 
			 rte_gp_lat, rte_gp_lon, i_space, ground_emission, 
			 ground_los, ground_refl_coeffs, f_index, scat_za_index, 
			 scat_aa_index, ext_mat_spt, abs_vec_spt, pnd_ppathLOS, 
			 ppath_step_agenda, atmosphere_dim, 
			 p_grid, lat_grid, lon_grid, z_field, r_geoid, z_ground, 
			 cloudbox_limits, record_ppathcloud, record_ppath, 
			 opt_prop_gas_agenda, spt_calc_agenda, 
			 scalar_gas_absorption_agenda, stokes_dim, t_field, 
			 vmr_field, rte_agenda, i_space_agenda, 
			 ground_refl_agenda, f_grid, 0, 0,pnd_field);
  

  mult(IboundaryLOScontri,TArrayLOS[TArrayLOS.nelem()-1],i_rte(0,joker));
  for (Index i = 0;i<stokes_dim;i++){assert(!isnan(IboundaryLOScontri[i]));}
  
  //Begin Main Loop
  for (photon_number=1; photon_number<=maxiter; photon_number++)
    {
      keepgoing=true; //flag indicating whether to continue tracing a photon path
      scattering_order=0;	//scattering order
      Q=identity;		//identity matrix
      pathI=0.0;
      boundarycontri=0.0;
      pathinc=0.0;
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
	 
      if (silent==0){cout<<"photon_number = "<<photon_number<<"\n";}
      while (keepgoing)
	{
	   if (scattering_order>0)
	    {
	      //We need to calculate a new propagation path. In the future, we may be 
	      //able to take some shortcuts here
	      Cloudbox_ppathCalc(ppathcloud,ppath_step,ppath_step_agenda,atmosphere_dim,
				 p_grid,lat_grid,lon_grid,z_field,r_geoid,z_ground,
				 cloudbox_limits, rte_pos,rte_los);
	      if (record_ppathcloud){ppathRecordMC(ppathcloud,"ppathcloud",
						   photon_number,scattering_order);}
			      
	      cum_l_stepCalc(cum_l_step,ppathcloud);
  
	      //Calculate array of transmittance matrices
	      TArrayCalc(TArray, ext_matArray, abs_vecArray, t_ppath, scat_za_grid, 
			 scat_aa_grid, ext_mat, abs_vec, rte_pressure, rte_temperature, 
			 rte_vmr_list, scat_za_index, scat_aa_index, ext_mat_spt, 
			 abs_vec_spt, pnd_ppath, ppathcloud, opt_prop_gas_agenda, 
			 spt_calc_agenda, scalar_gas_absorption_agenda, stokes_dim, 
			 p_grid, lat_grid, lon_grid, t_field, vmr_field, atmosphere_dim,
			 pnd_field);


	      /////////////////////////////////////////////////////////////////////
	      dist_to_boundary=cum_l_step[ppathcloud.np-1];
	 
	      //	      Iboundary=i_rte(0,joker);
	      Sample_ppathlength (pathlength,g,rng,ext_matArray,TArray,cum_l_step,2);
	    }
	  else
	    {
	      Sample_ppathlengthLOS (pathlength,g,rng,ext_matArray,TArray,cum_l_step,2);
	    }
	  assert(cum_l_step.nelem()==ppathcloud.np);
	  assert(TArray.nelem()==ppathcloud.np);
	  if (pathlength>dist_to_boundary)
	    //Then the path has left the cloud box
	    {
	      assert (scattering_order>0); //scattering/emission should be 
				           //forced in original line of sight
	      //Get incoming//////
	      montecarloGetIncoming(i_rte,rte_pos,rte_los,rte_gp_p,
			rte_gp_lat,rte_gp_lon,ppath,ppath_step,i_space,
			ground_emission,ground_los,ground_refl_coeffs,
			scat_za_grid,scat_aa_grid,ppath_step_agenda,
			rte_agenda,i_space_agenda,ground_refl_agenda,t_field,
			p_grid,lat_grid,lon_grid,z_field,r_geoid,
			z_ground,cloudbox_limits,ppathcloud,atmosphere_dim,
			f_grid,stokes_dim);
	      
	      if (record_ppath)
		{
		  ppathRecordMC(ppath,"ppath",photon_number,scattering_order);
		}
	      
	      f_index=0;//For some strange reason f_index is set to -1 in RteStandard
	      Iboundary=i_rte(0,joker);
	      ////////////////////
	      T=TArray[ppathcloud.np-1];
	      mult(boundarycontri,T,Iboundary);
	      mult(pathinc,Q,boundarycontri);
	      pathI = pathinc;
	      pathI/=g;//*=exp(K11*dist_to_boundary);
	      if(strat_sampling)
		{
		  p_from_below=p_from_belowCscaCalc(za_scat,
						    montecarlo_p_from_belowCsca,
						    pnd_vec,za_grid)/(K(0,0)-K_abs[0]);
		  
		  if(rte_los[0]>90)
		    {
		      pathI*=p_from_below;//Needs to be renormalised
		      I_from_below_sum+=pathI;
		      I_from_below_N+=1;
		      for(Index j=0; j<stokes_dim; j++)
			{
			  assert(!isnan(pathI[j]));
			  I_from_below_squaredsum[j] += pathI[j]*pathI[j];
			}
		    }
		  else
		    {
		      pathI*=(1-p_from_below);//Needs to be renormalised
		      I_from_above_sum+=pathI;
		      I_from_above_N+=1;
		      for(Index j=0; j<stokes_dim; j++)
			{
			  assert(!isnan(pathI[j]));
			  I_from_above_squaredsum[j] += pathI[j]*pathI[j];
			}
		    }
		      
		}
	      keepgoing=false; //stop here. New photon.
	    }
	  else
	    {
	       //we have another scattering/emission point
	      //Interpolate T, s_0, etc from ppath and Tarray
	      interpTArray(T, K_abs, temperature, K, rte_pos, rte_los, pnd_vec,
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
		  Numeric planck_value = planck( f_grid[f_index], temperature );
		  Vector emission=K_abs;
		  emission*=planck_value;
		  Vector emissioncontri(stokes_dim);
		  mult(emissioncontri,T,emission);
		  emissioncontri/=(g*(1-albedo));//yuck!
		  mult(pathinc,Q,emissioncontri);
		  pathI = pathinc;
		  if(strat_sampling)
		    {
		      pathI*=1-albedo;
		      I_emission_sum+=pathI;
		      I_emission_N+=1;
		      for(Index j=0; j<stokes_dim; j++)
			{
			  assert(!isnan(pathI[j]));
			  I_emission_squaredsum[j] += pathI[j]*pathI[j];
			}
		    }
		  keepgoing=false;
		}
	      else
		{
		  //Sample new line of sight.
		  if (los_sampling_method==3)
		    {
		      Sample_los_Z (new_rte_los,g_los_csc_theta,Z,rng,rte_los,
				    scat_data_raw,stokes_dim,f_index,f_grid,
				    scat_theta,scat_theta_gps,scat_theta_itws,
				    pnd_vec,K(0,0)-K_abs[0]);
		    }
		  else
		    {
		      Sample_los(new_rte_los,g_los_csc_theta,rng, los_sampling_method);
		      
		      //Calculate Phase matrix////////////////////////////////
		      Numeric aa_scat = (rte_los[1]>=0) ?-180+rte_los[1]:180+rte_los[1];
		      Numeric aa_inc= (new_rte_los[1]>=0) ?
			-180+new_rte_los[1]:180+new_rte_los[1];
		      
		      pha_mat_singleCalc(Z,180-rte_los[0],aa_scat,180-new_rte_los[0],
					 aa_inc,scat_data_raw,stokes_dim,f_index,
					 f_grid,scat_theta,scat_theta_gps,
					 scat_theta_itws,pnd_vec);
		    }
		  if(strat_sampling)
		    {
		      Z/=g*g_los_csc_theta;
		    }
		  else
		    {
		      Z/=g*g_los_csc_theta*albedo;
		    }
		  mult(q,T,Z);
		  mult(newQ,Q,q);
		  Q=newQ;
		  scattering_order+=1;
		  za_scat=180-rte_los[0];
		  rte_los=new_rte_los;
		  if (silent==0){cout <<"photon_number = "<<photon_number << 
				   ", scattering_order = " <<scattering_order <<"\n";}
		}

	    }
 
	}
      if (!strat_sampling)
      {
	Isum += pathI;
	if (record_histdata==1){histfile << pathI << "\n";}
	for(Index j=0; j<stokes_dim; j++)
	  {
	    assert(!isnan(pathI[j]));
	    Isquaredsum[j] += pathI[j]*pathI[j];
	  }
      }
      if (photon_number==500)
	{
	  cout <<"Estimated execution time for ScatteringMonteCarlo: " << 
	    (Numeric)(time(NULL)-start_time)*maxiter/500 <<" seconds.\n";
	}
    }
  if (strat_sampling)
    {
      Vector I_emission=I_emission_sum;
      I_emission/=I_emission_N;
      //cout<<I_emission<<"\n";
      Vector I_from_below=I_from_below_sum;
      I_from_below/=I_from_below_N;
      //cout<<I_from_below<<"\n";
      Vector I_from_above=I_from_above_sum;
      I_from_above/=I_from_above_N;
      //cout<<I_from_above<<"\n";
      for(Index j=0; j<stokes_dim; j++)	
	{
	  I[j]=I_emission[j]+I_from_below[j]+I_from_above[j];
	  i_montecarlo_error[j]=sqrt((I_from_below_squaredsum[j]/I_from_below_N-
				      I_from_below[j]*I_from_below[j])
				     /I_from_below_N+
				     (I_from_above_squaredsum[j]/I_from_above_N-
				      I_from_above[j]*I_from_above[j])
				     /I_from_above_N+
				     (I_emission_squaredsum[j]/I_emission_N-
				      I_emission[j]*I_emission[j])
				     /I_emission_N);
	}
    }
  else
    {
      I=Isum;
      I/=maxiter;
      for(Index j=0; j<stokes_dim; j++)	
	{
	  i_montecarlo_error[j]=sqrt((Isquaredsum[j]/maxiter-I[j]*I[j])/maxiter);
	}
    }
  I+=IboundaryLOScontri;
  i_rte(0,joker)=I;
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





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
#include <sstream>
#include "rte.h"
#include "lin_alg.h"
#include "auto_md.h"
#include "logic.h"
#include "physics_funcs.h"
#include "xml_io.h"
#include "montecarlo.h"
#include "rng.h"


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
        const Vector&         a_pos,
        const Vector&         a_los)
{
  // This function is a WSM but it is normally only called from RteCalc. 
  // For that reason, this function does not repeat input checks that are
  // performed in RteCalc, it only performs checks regarding the sensor 
  // position and LOS.

  //--- Check input -----------------------------------------------------------

  // Sensor position and LOS
  //
  chk_vector_length( "a_pos", a_pos, atmosphere_dim );
  chk_if_over_0( "sensor radius", a_pos[0] );
  if( atmosphere_dim < 3 )
    {
	ostringstream os;
	os << "cloudbox_ppath_calc only works for a 3D atmosphere";
	throw runtime_error( os.str() );
    }
  else
    {
      chk_if_in_range( "sensor latitude", a_pos[1], -90, 90 );
      chk_if_in_range( "sensor longitude", a_pos[2], -360, 360 );
      chk_vector_length( "a_los", a_los, 2 );
      chk_if_in_range( "sensor zenith angle", a_los[0], 0, 180 );
      chk_if_in_range( "sensor azimuth angle", a_los[1], -180, 180 );
    }
  
  //--- End: Check input ------------------------------------------------------


  // Some messages
  out2 << "  -------------------------------------\n";
  out2 << "  sensor radius          : " << a_pos[0]/1e3 << " km\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor longitude       : " << a_pos[2] << "\n";
  out2 << "  sensor zenith angle    : " << a_los[0] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor azimuth angle   : " << a_los[1] << "\n";
  
  
  
  // Initiate the partial Ppath structure. 

  //
  cloudbox_ppath_start_stepping( ppath_step, atmosphere_dim, p_grid, lat_grid, 
                        lon_grid, z_field, r_geoid, z_ground,a_pos, a_los );

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
			   Vector&               a_pos,
			   Vector&               a_los,
			   //Stuff needed by RteCalc
			   GridPos&              a_gp_p,
			   GridPos&              a_gp_lat,
			   GridPos&              a_gp_lon,
			   Matrix&               i_space,
			   Matrix&               ground_emission,
			   Matrix&               ground_los, 
			   Tensor4&              ground_refl_coeffs,
			   Matrix&               i_rte,
			   Vector&               scat_za_grid,
			   Vector&               scat_aa_grid,
			   Numeric&              a_pressure,
			   Numeric&              a_temperature,
			   Vector&               a_vmr_list,
			   //Other Stuff
			   Tensor3&              ext_mat,
			   Matrix&               abs_vec,
			   Index&                f_index,
			   Index&                scat_za_index,
			   Index&                scat_aa_index,
			  
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
			   const Agenda& opt_prop_part_agenda,
			   const Agenda& scalar_gas_absorption_agenda,
			   const Tensor4&   vmr_field,
			   //Other Stuff
			   const ArrayOfSingleScatteringData& scat_data_raw,
			   const Tensor4& pnd_field,
			    // Control Parameters:
			   const Index& maxiter,
			   const Index& rng_seed,
			   const Index& record_ppathcloud,
			   const Index& record_ppath)

{		
  //Internal Declarations
   Matrix identity(stokes_dim,stokes_dim,0.0);
  //Identity matrix
  for (Index i=0; i<stokes_dim; i++){identity(i,i)=1.0;}
  Matrix Q(stokes_dim,stokes_dim),T(stokes_dim,stokes_dim);
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim);
  bool keepgoing;
  Index scattering_order;
  Vector new_a_los(2);
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
  Index za_prop;
  Index aa_prop;
  Tensor5 pha_mat_spt(1,2,2,stokes_dim,stokes_dim);
  Tensor4 pha_mat(2,2,stokes_dim,stokes_dim);
  Index p_index;
  Index lat_index;
  Index lon_index;
   ArrayOfGridPos pathlength_gp(1);
  Vector K_abs(stokes_dim);
  Vector I(stokes_dim);
  i_montecarlo_error.resize(stokes_dim);
  Rng rng;

  //if rng_seed is < 0, keep time based seed, otherwise...
  if(rng_seed>=0){rng.seed(rng_seed);}
  
  Cloudbox_ppath_rteCalc(ppathLOS, ppath, ppath_step, a_pos, a_los, cum_l_stepLOS, 
		 TArrayLOS, ext_matArrayLOS, abs_vecArrayLOS,t_ppathLOS, scat_za_grid, 
			 scat_aa_grid, ext_mat, abs_vec, a_pressure, a_temperature, 
			 a_vmr_list, i_rte, a_gp_p, a_gp_lat, a_gp_lon, 
			 i_space, ground_emission, ground_los, ground_refl_coeffs, 
			 f_index, scat_za_index, scat_aa_index,
			 ppath_step_agenda, atmosphere_dim, p_grid, 
			 lat_grid, lon_grid, z_field, r_geoid, z_ground, 
			 cloudbox_limits, record_ppathcloud, record_ppath, 
			 opt_prop_gas_agenda, opt_prop_part_agenda, 
			 scalar_gas_absorption_agenda, stokes_dim, t_field, vmr_field, 
			 rte_agenda, i_space_agenda, ground_refl_agenda, f_grid, 0, 0);
  

  mult(IboundaryLOScontri,TArrayLOS[TArrayLOS.nelem()-1],i_rte(0,joker));
  
  //Begin Main Loop
  for (Index photon_number=1; photon_number<=maxiter; photon_number++)
    {
      keepgoing=true;      //flag indicating whether to continue tracing a photon path
      scattering_order=0;              //scattering order
      Q=identity;       //identity matrix
      Vector pathI(stokes_dim,0.0);
      Vector boundarycontri(stokes_dim,0.0);
      Vector pathinc(stokes_dim,0.0);
      //while the reversed traced photon path remains in the cloud box
      //
      TArray=TArrayLOS;
      ext_matArray=ext_matArrayLOS;
      abs_vecArray=abs_vecArrayLOS;
      ppathcloud=ppathLOS;
      cum_l_step=cum_l_stepLOS;
      t_ppath=t_ppathLOS;
      cout<<"photon_number = "<<photon_number<<"\n";
      while (keepgoing)
	{
	  if (scattering_order>0)
	    {
	      //We need to calculate a new propagation path. In the future, we may be 
	      //able to take some shortcuts here
	      Cloudbox_ppath_rteCalc(ppathcloud, ppath, ppath_step, a_pos, a_los, cum_l_step, 
				    TArray, ext_matArray, abs_vecArray,t_ppath, scat_za_grid, 
				     scat_aa_grid, ext_mat, abs_vec, a_pressure, 
				     a_temperature, a_vmr_list, i_rte, a_gp_p, 
				     a_gp_lat, a_gp_lon, i_space, ground_emission, 
				     ground_los, ground_refl_coeffs, f_index, scat_za_index, 
				     scat_aa_index, 
				     ppath_step_agenda, atmosphere_dim, p_grid, lat_grid, 
				     lon_grid, z_field, r_geoid, z_ground, cloudbox_limits, 
				     record_ppathcloud, record_ppath, opt_prop_gas_agenda, 
				     opt_prop_part_agenda, scalar_gas_absorption_agenda, 
				     stokes_dim, t_field, vmr_field, rte_agenda, 
				     i_space_agenda, ground_refl_agenda, f_grid, photon_number
				     , scattering_order);
	      Iboundary=i_rte(0,joker);
	    }
	  
	  Sample_ppathlength (pathlength,g,rng,ext_matArray);
	  //	  cout << "pathlength = " << pathlength << "\n";


	  assert(cum_l_step.nelem()==ppathcloud.np);
	  assert(TArray.nelem()==ppathcloud.np);
	  if (pathlength>cum_l_step[ppathcloud.np-1])
	    //Then the path has left the cloud box
	    {
	      if (scattering_order>0)
		{
		  T=TArray[ppathcloud.np-1];
		  mult(boundarycontri,T,Iboundary);
		  mult(pathinc,Q,boundarycontri);
		  pathI += pathinc;
		}
	      keepgoing=false;     
	      //stop here
	    }
	  else
	    {
	       //we have another scattering/emission point
	      //Interpolate T, s_0, etc from ppath and Tarray
	      interpTArray(T, K_abs, temperature, a_pos, a_los, pathlength_gp,TArray, 
			   ext_matArray,abs_vecArray, t_ppath, cum_l_step,pathlength, 
			   stokes_dim, ppathcloud);
	      //Calculate emission
	      Numeric planck_value = planck( f_grid[f_index], temperature );
	      Vector emission=K_abs;
	      emission*=planck_value;
	      Vector emissioncontri(stokes_dim);
	      mult(emissioncontri,T,emission);
	      emissioncontri/=g;
	      if (scattering_order>0)
		{
		  mult(boundarycontri,TArray[ppathcloud.np-1],Iboundary);
		  emissioncontri+=boundarycontri;  
		}
	      mult(pathinc,Q,emissioncontri);
	      pathI += pathinc;
	      Sample_los(new_a_los,rng);
	      //Calculate Phase matrix////////////////////////////////
	      pha_mat_za_grid[0]=180-a_los[0];
	      pha_mat_za_grid[1]=180-new_a_los[0];
	      pha_mat_aa_grid[0]=180+a_los[1];
	      pha_mat_aa_grid[1]=180+new_a_los[1];
	      za_prop=0;
	      aa_prop=0;
	      pha_mat_sptFromData(pha_mat_spt,
				  scat_data_raw, pha_mat_za_grid, pha_mat_aa_grid, 
				  za_prop, aa_prop, f_index, f_grid);
	      //There should probably be some interpolation here for now just use the 
	      //low gridpoint
	      //  cout << "a_pos = "<<a_pos<<"\n";
	      
	      p_index=ppathcloud.gp_p[pathlength_gp[0].idx].idx;
	      lat_index=ppathcloud.gp_lat[pathlength_gp[0].idx].idx;
	      lon_index=ppathcloud.gp_lon[pathlength_gp[0].idx].idx;

	      pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
			  atmosphere_dim, p_index, lat_index, 
			  lon_index);
	      Z=pha_mat(1,1,joker,joker);
	      Z*=4*PI/g/(1-cos(DEG2RAD*a_los[0]));
	      mult(q,T,Z);
	      mult(newQ,Q,q);
	      Q=newQ;
	      scattering_order+=1;
	      a_los=new_a_los;
	      cout <<"photon_number = "<<photon_number << 
		", scattering_order = " <<scattering_order <<"\n";
	      if (Q(0,0)<1e-6){ keepgoing=false;}
	    }
	  //	  cout<<"pathI = "<<pathI<<"\n";
	  
	}
      Isum += pathI;
      for(Index j=0; j<stokes_dim; j++)
	{
	  Isquaredsum[j] += pathI[j]*pathI[j];
	}
    }
  I=Isum;
  I/=maxiter;
  I+=IboundaryLOScontri;
  for(Index j=0; j<stokes_dim; j++)	
    {
      i_montecarlo_error[j]=sqrt((Isquaredsum[j]/maxiter-I[j]*I[j])/maxiter);
    }
  i_rte(0,joker)=I;
}		




//! shift_a_pos
/*! 
   shifts a_pos and a_los to the end of ppath.  3D only!

   \author Cory Davis
   \date   2003-07-19
*/


void shift_a_pos(
	Vector&         a_pos,
      Vector&         a_los,
	const Ppath&    ppath)
{
a_pos = ppath.pos(ppath.np-1,Range(0,3));
a_los = ppath.los(ppath.np-1,joker);
}  




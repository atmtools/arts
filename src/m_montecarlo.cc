/* Copyright (C) 2003-2012 Cory Davis <cory@met.ed.ac.uk>
                            
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
#include "auto_md.h"
#include "ppath.h"
#include "matpackI.h"
#include "special_interp.h"
#include "check_input.h"
#include <stdexcept>
#include <cmath>
#include "rte.h"
#include "lin_alg.h"
#include "logic.h"
#include "physics_funcs.h"
#include "xml_io.h"
#include "montecarlo.h"
#include "rng.h"
#include <ctime>
#include <fstream>
#include "mc_interp.h"
#include "math_funcs.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;
extern const Numeric BOLTZMAN_CONST;
extern const Numeric SPEED_OF_LIGHT;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetGaussian(MCAntenna& mc_antenna,
                           //keyword arguments
                           const Numeric& za_sigma,
                           const Numeric& aa_sigma,
                           const Verbosity&)
{
  mc_antenna.set_gaussian(za_sigma,aa_sigma);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetGaussianByFWHM(MCAntenna& mc_antenna,
                                 //keyword arguments
                                 const Numeric& za_fwhm,
                                 const Numeric& aa_fwhm,
                                 const Verbosity&)
{
  mc_antenna.set_gaussian_fwhm(za_fwhm,aa_fwhm);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetPencilBeam(MCAntenna& mc_antenna,
                             const Verbosity&)
{
  mc_antenna.set_pencil_beam();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MCGeneral(Workspace&            ws,
               Vector&               y,
               Index&                mc_iteration_count,
               Vector&               mc_error,
               Tensor3&              mc_points,
               ArrayOfIndex&         mc_source_domain,
               ArrayOfIndex&         mc_scat_order,
               const MCAntenna&      mc_antenna,
               const Vector&         f_grid,
               const Index&          f_index,
               const Matrix&         sensor_pos,
               const Matrix&         sensor_los,
               const Index&          stokes_dim,
               const Index&          atmosphere_dim,
               const Agenda&         ppath_step_agenda,
               const Numeric&        ppath_lmax,
               const Numeric&        ppath_lraytrace,
               const Agenda&         iy_space_agenda,
               const Agenda&         surface_rtprop_agenda,
               const Agenda&         propmat_clearsky_agenda, 
               const Vector&         p_grid,
               const Vector&         lat_grid, 
               const Vector&         lon_grid, 
               const Tensor3&        z_field, 
               const Vector&         refellipsoid, 
               const Matrix&         z_surface,
               const Tensor3&        t_field, 
               const Tensor4&        vmr_field, 
               const Index&          cloudbox_on,
               const ArrayOfIndex&   cloudbox_limits, 
               const Tensor4&        pnd_field,
               const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
               const Index&          atmfields_checked,
               const Index&          atmgeom_checked,
               const Index&          cloudbox_checked,
               const String&         iy_unit,
               const Index&          mc_seed,
               const Numeric&        std_err,
               const Index&          max_time,
               const Index&          max_iter,
               const Index&          min_iter,
               const Numeric&        taustep_limit,
               const Index&          l_mc_scat_order,
               const Verbosity&      verbosity)
{
  // Checks of input
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );
  if( !cloudbox_on )
    throw runtime_error( "The cloudbox  must be activated (cloudbox_on=1)." );

  if( min_iter < 100 )
    throw runtime_error( "*mc_min_iter* must be >= 100." );

  if( max_time < 0  &&  max_iter < 0  &&  std_err < 0 )
    throw runtime_error( "At least one of std_err, max_time, and max_iter "
                         "must be positive." );

  if( l_mc_scat_order <= 0 )
    throw runtime_error( "*l_max_scat_order* must be > 0." );

  if( f_index < 0 )
    throw runtime_error( "The option of f_index < 0 is not handled by this "
                         "method." );
  if( f_index >= f_grid.nelem() )
    throw runtime_error( "*f_index* is outside the range of *f_grid*." );

  if( atmosphere_dim != 3 )
    throw runtime_error( "Only 3D atmospheres are handled. " );

  if( sensor_pos.ncols() != 3 )
    {
      ostringstream os;
      os << "Expected number of columns in sensor_pos: 3.\n";
      os << "Found: " << sensor_pos.ncols();
      throw runtime_error(os.str());
    }

  if (! (sensor_los.ncols() == 2))
    {
      ostringstream os;
      os << "Expected number of columns in sensor_los: 2.\n";
      os << "Found: " << sensor_los.ncols();
      throw runtime_error(os.str());
    }

  for( Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss ++ )
    {
      for( Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se ++ )
        {
          if( scat_data_mono[i_ss][i_ss].f_grid.nelem() > 1 )
            {
              ostringstream os;
              os << "Input scattering data shall be monochromatic, but at least\n"
                 << "one element of scat_data_mono contains data for several\n"
                 << "frequencies.\n"
                 << "Found incorrect element is [" << i_ss << "][" << i_se 
                 << "] (zero based indexing).";
              throw runtime_error( os.str() );
            }   
          if( abs( ( scat_data_mono[i_ss][i_ss].f_grid[0] - f_grid[f_index] ) /
                   f_grid[f_index] ) > 1e-4 )
            {
              ostringstream os;
              os << "Incorrect frequency in at least one element of scat_data_mono.\n" 
                 << "Scattering data have frequency : " << 
                           scat_data_mono[i_ss][i_ss].f_grid[0]/1e9 << " GHz\n"
                 << "Radiative transfer frequency is: " << 
                           f_grid[f_index]/1e9 << " GHz\n"
                 << "Found incorrect element is [" << i_ss << "][" << i_se 
                 << "] (zero based indexing).\n"
                 << "Maximum allowed relative frequency is 1e-4.";
              throw runtime_error( os.str() );
            }   
        }   
    }


  Ppath  ppath_step;
  Rng    rng;                      //Random Number generator
  time_t start_time=time(NULL);
  Index  N_se = pnd_field.nbooks();//Number of scattering elements
  Vector pnd_vec(N_se); //Vector of particle number densities used at each point
  Vector Z11maxvector;//Vector holding the maximum phase function for each 
  bool  anyptype30 = is_anyptype30(scat_data_mono);
  if (anyptype30)
    { findZ11max(Z11maxvector,scat_data_mono); }
  rng.seed(mc_seed, verbosity);
  Numeric g,temperature,albedo,g_los_csc_theta;
  Matrix  A(stokes_dim,stokes_dim), Q(stokes_dim,stokes_dim);
  Matrix  evol_op(stokes_dim,stokes_dim), ext_mat_mono(stokes_dim,stokes_dim);
  Matrix  q(stokes_dim,stokes_dim), newQ(stokes_dim,stokes_dim);
  Matrix  Z(stokes_dim,stokes_dim);
  q     = 0.0; 
  newQ  = 0.0;
  Vector vector1(stokes_dim), abs_vec_mono(stokes_dim), I_i(stokes_dim);
  Vector Isum(stokes_dim), Isquaredsum(stokes_dim);
  Index termination_flag = 0;
  const Numeric f_mono = f_grid[f_index];

  CREATE_OUT0;

  y.resize(stokes_dim);
  y = 0;

  mc_iteration_count = 0;
  mc_error.resize(stokes_dim);
  mc_points.resize( p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem() );
  mc_points = 0;
  mc_scat_order.resize( l_mc_scat_order );
  mc_scat_order = 0;
  mc_source_domain.resize( 4 );
  mc_source_domain = 0;


  //local versions of workspace
  Matrix  local_iy(1,stokes_dim), local_surface_emission(1,stokes_dim);
  Matrix  local_surface_los;
  Tensor4 local_surface_rmatrix;
  Vector  local_rte_pos(2);
  Vector  local_rte_los(2);
  Vector  new_rte_los(2);
  Index   np;
  Isum=0.0;Isquaredsum=0.0;
  Numeric std_err_i;
  bool    convert_to_rjbt = false;
  if ( iy_unit == "RJBT" )
    { 
      std_err_i = f_mono*f_mono*2*BOLTZMAN_CONST/
                                          SPEED_OF_LIGHT/SPEED_OF_LIGHT*std_err;
      convert_to_rjbt = true;
    }
  else if ( iy_unit == "1" )
    { std_err_i = std_err; }
  else
    {
      ostringstream os;
      os << "Invalid value for *iy_unit*:" << iy_unit <<".\n" 
         << "This method allows only the options \"RJBT\" and \"1\".";
      throw runtime_error( os.str() );
    }
      

  //Begin Main Loop
  //
  bool keepgoing, oksampling;
  Index nfails = 0;
  //
  while( true )
    {
      // Complete content of while inside try/catch to handle occasional
      // failuers in the ppath calculations
      try{
        
        bool inside_cloud;

        mc_iteration_count += 1;
        Index scattering_order = 0;

        keepgoing = true;    // indicating whether to continue tracing a photon
        oksampling = true;   // gets false if g becomes zero

        //Sample a FOV direction
        mc_antenna.draw_los(local_rte_los,rng,sensor_los(0,joker));
        id_mat(Q);
        local_rte_pos=sensor_pos(0,joker);
        I_i=0.0;
      
        while( keepgoing )
          {
            mcPathTraceGeneral( ws, evol_op, abs_vec_mono, temperature, 
                                ext_mat_mono, rng, local_rte_pos, local_rte_los, 
                                pnd_vec, g, ppath_step, termination_flag, 
                                inside_cloud, ppath_step_agenda,
                                ppath_lmax, ppath_lraytrace, taustep_limit,
                                propmat_clearsky_agenda, stokes_dim, f_mono, 
                                p_grid, lat_grid, lon_grid, z_field, refellipsoid,
                                z_surface, t_field, vmr_field, 
                                cloudbox_limits, pnd_field, scat_data_mono, 
                                verbosity ); 
             

            // GH 2011-09-08: if the lowest layer has large
            // extent and a thick cloud, g may be 0 due to
            // underflow, but then I_i should be 0 as well.
            // Don't turn it into nan for no reason.
            // If reaching underflow, no point in going on;
            // hence new photon.
            // GH 2011-09-14: moved this check to outside the different
            // scenarios, as this goes wrong regardless of the scenario.
            if( g == 0 )
              {
                keepgoing  = false;
                oksampling = false;              
                mc_iteration_count -= 1;
                out0 << "WARNING: A rejected path sampling (g=0)!\n(if this"
                     << "happens repeatedly, try to decrease *ppath_lmax*)";
              }
            else if( termination_flag == 1 )
              {
                iy_space_agendaExecute( ws, local_iy, Vector(1,f_mono), 
                                        local_rte_pos, local_rte_los, 
                                        iy_space_agenda );
                mult( vector1, evol_op, local_iy(0,joker) );
                mult( I_i, Q, vector1 );
                I_i /= g;
                keepgoing=false; //stop here. New photon.
                mc_source_domain[0] += 1;
              }
            else if( termination_flag == 2 )
              {
                //Calculate surface properties
                surface_rtprop_agendaExecute( ws, local_surface_emission, 
                                              local_surface_los, 
                                              local_surface_rmatrix, 
                                              Vector(1,f_mono),
                                              local_rte_pos, local_rte_los, 
                                              surface_rtprop_agenda );
                
                //if( local_surface_los.nrows() > 1 )
                // throw runtime_error( 
                //                "The method handles only specular reflections." );

                //deal with blackbody case
                if( local_surface_los.empty() )
                  {
                    mult( vector1, evol_op, local_surface_emission(0,joker) );
                    mult( I_i, Q, vector1);
                    I_i /= g;
                    keepgoing = false;
                    mc_source_domain[1] += 1;
                  }
                else
                  //decide between reflection and emission
                  {
                    const Numeric rnd = rng.draw();

                    Numeric R11 = 0;
                    for( Index i=0; i<local_surface_rmatrix.nbooks(); i++ )
                      { R11 += local_surface_rmatrix(i,0,0,0); }
                    
                    if( rnd > R11 )
                      {
                        //then we have emission
                        mult( vector1, evol_op, local_surface_emission(0,joker) );
                        mult( I_i, Q, vector1 );
                        I_i /= g*(1-R11);
                        keepgoing = false;
                        mc_source_domain[1] += 1;
                      }
                    else
                      {
                        //we have reflection
                        // determine which reflection los to use
                        Index i=0;
                        Numeric rsum = local_surface_rmatrix(i,0,0,0);
                        while( rsum < rnd )
                          {
                            i++;
                            rsum += local_surface_rmatrix(i,0,0,0);
                          }
                        
                        local_rte_los = local_surface_los( i, joker );
                        
                        mult( q, evol_op, local_surface_rmatrix(i,0,joker,joker) );
                        mult( newQ, Q, q );
                        Q  = newQ;
                        Q /= g * local_surface_rmatrix(i,0,0,0);
                      }
                  }
              }
            else if (inside_cloud)
              {
                //we have another scattering/emission point 
                //Estimate single scattering albedo
                albedo = 1 - abs_vec_mono[0]/ext_mat_mono(0,0);

                //determine whether photon is emitted or scattered
                if( rng.draw() > albedo )
                  {
                    //Calculate emission
                    Numeric planck_value = planck( f_mono, temperature );
                    Vector emission = abs_vec_mono;
                    emission *= planck_value;
                    Vector emissioncontri(stokes_dim);
                    mult( emissioncontri, evol_op, emission );
                    emissioncontri /= (g*(1-albedo)); //yuck!
                    mult( I_i, Q, emissioncontri );
                    keepgoing = false;
                    mc_source_domain[3] += 1;                  
                  }
                else
                  {
                    //we have a scattering event
                    Sample_los( new_rte_los, g_los_csc_theta, Z, rng, 
                                local_rte_los, scat_data_mono, stokes_dim,
                                pnd_vec, anyptype30, Z11maxvector, 
                                ext_mat_mono(0,0)-abs_vec_mono[0], temperature,
                                verbosity );
                                             
                    Z /= g * g_los_csc_theta * albedo;
                    
                    mult( q, evol_op, Z );
                    mult( newQ, Q, q );
                    Q = newQ;
                    scattering_order += 1;
                    local_rte_los = new_rte_los;
                  }
              }
            else
              {
                //Must be clear sky emission point
                //Calculate emission
                Numeric planck_value = planck( f_mono, temperature );
                Vector emission = abs_vec_mono;
                emission *= planck_value;
                Vector emissioncontri(stokes_dim);
                mult( emissioncontri, evol_op, emission );
                emissioncontri /= g;
                mult( I_i, Q, emissioncontri );
                keepgoing = false;
                mc_source_domain[2] += 1;                  
              }
          }  // keepgoing

        if( oksampling )
          {
            // Set spome of the bookkeeping variables
            np = ppath_step.np;
            mc_points(ppath_step.gp_p[np-1].idx,ppath_step.gp_lat[np-1].idx,
                                                ppath_step.gp_lon[np-1].idx) += 1;
            if( scattering_order < l_mc_scat_order )
              { mc_scat_order[scattering_order] += 1; }

            Isum += I_i;

            for( Index j=0; j<stokes_dim; j++ )
              {
                assert( !isnan(I_i[j]) );
                Isquaredsum[j] += I_i[j]*I_i[j];
              }
            y  = Isum;
            y /= (Numeric)mc_iteration_count;
            for(Index j=0; j<stokes_dim; j++) 
              {
                mc_error[j] = sqrt( ( Isquaredsum[j] / 
                                      (Numeric)mc_iteration_count - y[j]*y[j] ) / 
                                      (Numeric)mc_iteration_count );
              }
            if( std_err > 0  &&  mc_iteration_count >= min_iter && 
                mc_error[0] < std_err_i )
              { break; }
            if( max_time > 0  &&  (Index)(time(NULL)-start_time) >= max_time )
              { break; }
            if( max_iter > 0  &&  mc_iteration_count >= max_iter )
              { break; }
          }
      }  // Try

      catch (runtime_error e)
        {
          mc_iteration_count += 1;
          nfails += 1;
          out0 << "WARNING: A MC path sampling failed! Error was:\n";
          cout << e.what () << endl;
          if( nfails >= 5 )
            {
              throw runtime_error(
              "The MC path sampling has failed five times. A few failures "
              "should be OK, but this number is suspiciously high and the "
              "reason to these failures should be tracked down." );
            }
        }
    }  // while

  if( convert_to_rjbt )
    {
      for(Index j=0; j<stokes_dim; j++) 
        {
          y[j]        = invrayjean( y[j],        f_mono );
          mc_error[j] = invrayjean( mc_error[j], f_mono );
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void MCSetSeedFromTime(Index& mc_seed,
                       const  Verbosity&)
{
  mc_seed=(Index)time(NULL);
}



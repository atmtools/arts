/* Copyright (C) 2006-2012 Claudia Emde <claudia.emde@dlr.de>
                      
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
   USA. 
*/

/**
 * \file   disort.cc
 * \author Claudia Emde <claudia.emde@dlr.de>
 * \date   Tue Feb  7 10:08:28 2006
 * 
 * \brief  This file contains functions related to the DISORT interface.   
 
**/

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "agenda_class.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "disort_DISORT.h"
#include "interpolation.h"
#include "logic.h"
#include "messages.h"
#include "rte.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric BOLTZMAN_CONST;

//! dtauc_ssalbCalc
/*!
  Calculates layer averaged cloud optical depth (dtauc) and 
  single scattering albedo (ssalb). These variables are required as
  input for the DISORT subroutine

  \param ws                    Current workspace
  \param dtauc                 optical depths for all layers
  \param ssalb                 single scattering albedos for all layers
  \param opt_prop_part_agenda  use arts -d for docu
  \param propmat_clearsky_agenda use arts -d
  \param spt_calc_agenda       use arts -d 
  \param pnd_field             use arts -d 
  \param t_field               use arts -d 
  \param z_field               use arts -d 
  \param p_grid                use arts -d 
  \param vmr_field             use arts -d 
  \param f_index               use arts -d 
  
  \author Claudia Emde, Jana Mendrok
  \date   2006-02-10
*/
void dtauc_ssalbCalc(Workspace& ws,
                     VectorView dtauc,
                     VectorView ssalb,
                     const Agenda& propmat_clearsky_agenda,
                     const Agenda& spt_calc_agenda,
                     const Agenda& opt_prop_part_agenda,
                     ConstTensor4View pnd_field,
                     ConstTensor3View t_field,
                     ConstTensor3View z_field, 
                     ConstTensor4View vmr_field,
                     ConstVectorView p_grid,
                     const ArrayOfIndex& cloudbox_limits,
                     ConstVectorView f_mono
                    )
{
  // Initialization
  dtauc=0.;
  ssalb=0.;
  
  const Index N_se = pnd_field.nbooks();
  // In DISORT the "cloudbox" must cover the whole atmosphere
  //const Index Np_cloud = pnd_field.npages();
  const Index Np = p_grid.nelem();

  //assert( dtauc.nelem() == Np_cloud-1);
  //assert( ssalb.nelem() == Np_cloud-1);
  assert( dtauc.nelem() == Np-1);
  assert( ssalb.nelem() == Np-1);

  const Index stokes_dim = 1; 

  // Local variables to be used in agendas
  Matrix abs_vec_spt_local(N_se, stokes_dim, 0.);
  Tensor3 ext_mat_spt_local(N_se, stokes_dim, stokes_dim, 0.);
  Matrix abs_vec_local;
  Tensor3 ext_mat_local;
  Numeric rtp_temperature_local; 
  Numeric rtp_pressure_local;
  Tensor4 propmat_clearsky_local;
  //Vector ext_vector(Np_cloud); 
  //Vector abs_vector(Np_cloud); 
  Vector ext_vector(Np, 0.); 
  Vector abs_vector(Np, 0.); 
  Vector rtp_vmr_local(vmr_field.nbooks());

  // Calculate ext_mat, abs_vec and sca_vec for all pressure points in cloudbox 
  for(Index scat_p_index_local = cloudbox_limits[0];
            scat_p_index_local < cloudbox_limits[1]+1; 
            scat_p_index_local ++)
    {
      rtp_temperature_local = t_field(scat_p_index_local, 0, 0);
     
      //Calculate optical properties for all individual scattering elements:
      spt_calc_agendaExecute(ws,
                             ext_mat_spt_local, 
                             abs_vec_spt_local,
                             scat_p_index_local, 0, 0, //position
                             rtp_temperature_local,
                             0, 0, // angles, only needed for za=0
                             spt_calc_agenda);

      opt_prop_part_agendaExecute(ws,
                                  ext_mat_local, abs_vec_local, 
                                  ext_mat_spt_local, 
                                  abs_vec_spt_local,
                                  scat_p_index_local, 0, 0, 
                                  opt_prop_part_agenda);

      ext_vector[scat_p_index_local] = ext_mat_local(0,0,0);
      abs_vector[scat_p_index_local] = abs_vec_local(0,0);
    }

  const Vector  rtp_temperature_nlte_local_dummy(0);

  // Calculate layer averaged single scattering albedo and layer optical depth
  propmat_clearsky_local = 0.;
  for (Index i = 0; i < Np-1; i++)
    {
      Numeric ext_part = 0.;
      Numeric abs_part = 0.;
 
      ext_part=.5*(ext_vector[i]+ext_vector[i+1]);
      abs_part=.5*(abs_vector[i]+abs_vector[i+1]);

      rtp_pressure_local = 0.5 * (p_grid[i] + p_grid[i+1]);
      rtp_temperature_local = 0.5 * (t_field(i,0,0) + t_field(i+1,0,0));
     
      // Average vmrs
      for (Index j = 0; j < vmr_field.nbooks(); j++)
        rtp_vmr_local[j] = 0.5 * (vmr_field(j, i, 0, 0) +
                                  vmr_field(j, i+1, 0, 0));
   
      const Vector rtp_mag_dummy(3,0);
      const Vector ppath_los_dummy;

      //FIXME: do this right?
      Tensor3 nlte_dummy;
      // This is right since there should be only clearsky partials
      ArrayOfTensor3 partial_dummy;
      // This is right since there should be only clearsky partials
      ArrayOfMatrix partial_source_dummy,partial_nlte_dummy;
      propmat_clearsky_agendaExecute(ws,
                                     propmat_clearsky_local,
                                     nlte_dummy,
                                     partial_dummy,
                                     partial_source_dummy,
                                     partial_nlte_dummy,
                                     ArrayOfRetrievalQuantity(0),
                                     f_mono,  // monochromatic calculation
                                     rtp_mag_dummy,ppath_los_dummy,
                                     rtp_pressure_local, 
                                     rtp_temperature_local, 
                                     rtp_temperature_nlte_local_dummy,
                                     rtp_vmr_local,
                                     propmat_clearsky_agenda);  

      //Assuming non-polarized light and only one frequency
      Numeric abs_gas = propmat_clearsky_local(joker,0,0,0).sum();

      if (ext_part!=0)
        ssalb[Np-2-i]=(ext_part-abs_part) / (ext_part+abs_gas);
     
      dtauc[Np-2-i]=(ext_part+abs_gas)*
        (z_field(i+1, 0, 0)-z_field(i, 0, 0));
    }  
}

//! phase_functionCalc
/*!
  Calculates layer averaged normalized phase functions from 
  the phase matrix in SingleScatteringData. The scattering angle 
  grid is taken from the data. 
  It is required that all scattering elements are given on the same 
  scattering angle grid (FIXME: Include angle interpolation)

  \param phase_function normalized phase function
  \param scat_data_mono as the WSV
  \param pnd_field      as the WSV
  
  \author Claudia Emde, Jana Mendrok
  \date   2006-02-10
*/
void phase_functionCalc(//Output
                        MatrixView phase_function,
                        //Input
                        const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                        ConstTensor4View pnd_field,
                        const ArrayOfIndex& cloudbox_limits)
{
/*
FIXME: dtauc_ssalbCalc applies spt_calc_agenda and opt_prop_part_agenda. Apply
 same/similar agendas here.
 It is inconsistent to consider T-dependence in dtauc_ssalbCalc, but neglect it
 here. It might be okayish, though, since T-dep of pfct might be lower than
 T-dep of ext/abs/scat coeffs.
*/

  // Initialization
  phase_function=0.;
  const Index nlyr = phase_function.nrows();

  const Index Np_cloud = pnd_field.npages();
  Matrix phase_function_level(Np_cloud, 
                              scat_data_mono[0][0].za_grid.nelem(), 0.);
  
  Vector sca_coeff_level(Np_cloud, 0.);

  //Loop over pressure levels
  for (Index i_p = 0; i_p < Np_cloud; i_p++)
    {
      // Calculate ensemble averaged scattering coefficient
      Numeric sca_coeff=0.;
      Index i_se_flat=0;
      //Numeric intP=0.;

      for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
        {
          for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
            {
              sca_coeff +=  pnd_field(i_se_flat, i_p, 0, 0) *
                (scat_data_mono[i_ss][i_se].ext_mat_data(0, 0, 0, 0, 0)-
                 scat_data_mono[i_ss][i_se].abs_vec_data(0, 0, 0, 0, 0));
              i_se_flat++;
            }
        }
      sca_coeff_level[i_p] = sca_coeff;

      // Bulk scattering function
      // (conversion to phase function only done when doing layer averaging.
      // this because averaging needs to be on scat coeff weighted phase
      // function aka bulk scattering function)
      if (sca_coeff != 0)
        {
          // Loop over scattering angles
          for (Index i_t = 0; i_t < scat_data_mono[0][0].za_grid.nelem(); i_t++)
            {
              i_se_flat=0;
              for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
                {
                  for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
                    {
                      phase_function_level(i_p, i_t) += 
                        pnd_field(i_se_flat, i_p, 0, 0) *
                        scat_data_mono[i_ss][i_se].pha_mat_data(0, 0, i_t,
                                                                0, 0, 0, 0);
                      i_se_flat++;
                    }
                }
            }
/*
              if( i_t>0 )
                  intP += PI *
                    (phase_function_level(i_p, i_t) +
                     phase_function_level(i_p, i_t-1)) *
                    abs(cos(scat_data_mono[0][0].za_grid[i_t]*PI/180.)-
                        cos(scat_data_mono[0][0].za_grid[i_t-1]*PI/180.));
            }

          cout << "at lev_cloud #" << i_p << ": "
               << "  total scatcoef=" << sca_coeff
               << ", integrated PFCT=" << intP << "\n";
*/
        }
    }


  // Calculate average phase function for the layers:
  // Average bulk scattering function and rescale (normalize) to phase function
  // with layer averaged scat coeff
  for (Index i_l = 0; i_l < Np_cloud-1; i_l++)
    {
      Index lyr_id = nlyr-1-i_l-cloudbox_limits[0];
      if ( phase_function_level(i_l, 0) !=0 )
        if( phase_function_level(i_l+1, 0) !=0 )
        {
          for (Index i_t=0; i_t < phase_function_level.ncols(); i_t++)
              phase_function(lyr_id, i_t) = 4*PI *
                ( phase_function_level(i_l, i_t) + 
                  phase_function_level(i_l+1, i_t) ) /
                ( sca_coeff_level[i_l] + sca_coeff_level[i_l+1] );
        }
        else
        {
          for (Index i_t=0; i_t < phase_function_level.ncols(); i_t++)
              phase_function(lyr_id, i_t) =
                phase_function_level(i_l, i_t) * 4*PI / sca_coeff_level[i_l];
        }
      else if ( phase_function_level(i_l+1, 0) !=0 )
      {
        for (Index i_t=0; i_t < phase_function_level.ncols(); i_t++)
            phase_function(lyr_id, i_t) =
              phase_function_level(i_l+1, i_t) * 4*PI / sca_coeff_level[i_l+1];
      }
    }
  
}

//! pmomCalc
/*!
  Calculates Legendre polynomials of phase functions for each layer. 
  The Legendre polynomial are required as input for DISORT. 

  \param pmom Legendre polynomial of phase functions
  \param phase_function Normalized phase function
  \param scat_angle_grid Scattering angle grid corresponding to phase 
  functions
  \param n_legendre Number of Legendre polynomials to be calculated
  
  \author Claudia Emde, Jana Mendrok
  \date   2006-02-10
*/
void pmomCalc(//Output
              MatrixView pmom,
              //Input
              ConstMatrixView phase_function, 
              ConstVectorView scat_angle_grid,
              const Index n_legendre,
              const Verbosity& verbosity)
{
  // Initialization
  pmom=0.;

  Numeric pint; //integrated phase function
  Numeric p0_1, p0_2, p1_1, p1_2, p2_1, p2_2;
  
  Vector za_grid(181);
  Vector u(181);

  for (Index i = 0; i< 181; i++)
    za_grid[i] = double(i);
  
  ArrayOfGridPos gp(181);
  gridpos(gp, scat_angle_grid, za_grid); 
  
  Matrix itw(gp.nelem(),2);    
  interpweights(itw,gp);
  
  Matrix phase_int(phase_function.nrows(),181);
  for  (Index i_l=0; i_l < phase_function.nrows(); i_l++)
    interp(phase_int(i_l, joker), itw, phase_function(i_l, joker), gp);
    
      
  for (Index i = 0; i<za_grid.nelem(); i++)
    u[i] = cos(za_grid[i] *PI/180.);
  
  for (Index i_l=0; i_l < phase_function.nrows(); i_l++)
    {
      pint = 0.;
      // Check if phase function is normalized
      for (Index i = 0; i<za_grid.nelem()-1; i++)            
        pint+=0.5*(phase_int(i_l, i)+phase_int(i_l, i+1))*
          abs(u[i+1] - u[i]);
      
      if (pint != 0){
        if (abs(2.-pint) > 0.2)
        {
          ostringstream os;
          os << "Phase function normalization deviates from expected value by\n"
             << "more than 20%. Something is wrong with your scattering data.\n"
             << "Check!\n";
          throw runtime_error( os.str() );
        }
        if (abs(2.-pint) > 1e-2)
        {
          CREATE_OUT2;
          out2 << "Warning: The phase function is not normalized to 2\n"
               << "The value is:" << pint << "\n";
        }
        
        //anyway, rescale phase_int to norm 2
        phase_int(i_l, joker) *= 2./pint;
       
        pmom(i_l, joker)= 0.; 

        for (Index i = 0; i<za_grid.nelem()-1; i++) 
          {
            p0_1=1.;
            p0_2=1.;
            
            pmom(i_l,0)=1.;

            //pmom(i_l,0)+=0.5*0.5*(phase_int(i_l, i)+phase_int(i_l, i+1))
            //*abs(u[i+1]-u[i]); 
            
            p1_1=u[i];
            p1_2=u[i+1];
            
            pmom(i_l,1)+=0.5*0.5*
              (p1_1*phase_int(i_l, i)+
               p1_2*phase_int(i_l, i+1))
              *abs(u[i+1]-u[i]);
            
            for (Index l=2; l<n_legendre; l++)
              {
              p2_1=(2*(double)l-1)/(double)l*u[i]*p1_1-((double)l-1)/
                (double)l*p0_1; 
              p2_2=(2*(double)l-1)/(double)l*u[i+1]*p1_2-((double)l-1)/
                (double)l*p0_2;
              
              pmom(i_l, l)+=0.5*0.5*
                (p2_1*phase_int(i_l, i)+
                 p2_2*phase_int(i_l, i+1))
                *abs(u[i+1]-u[i]);
              
              p0_1=p1_1;
              p0_2=p1_2;
              p1_1=p2_1;
              p1_2=p2_2;
              }
          }
        // cout << "pmom : " << pmom(i_l, joker) << endl;
        
      }
    }
}

#ifdef ENABLE_DISORT
 /* Workspace method: Doxygen documentation will be auto-generated */
void get_cb_inc_field(Workspace&      ws,
                      Matrix&         cb_inc_field,
                      const Agenda&   iy_main_agenda,
                      const Tensor3&  z_field,
                      const Tensor3&  t_field,
                      const Tensor4&  vmr_field,
                      const ArrayOfIndex&   cloudbox_limits,
                      const Vector&   f_grid,
                      const Vector&   scat_za_grid,
                      const Index&    nstreams
                     )
{
  // iy_unit hard.coded to "1" here
  const String iy_unit = "1";
  
  Matrix iy;

  //Define the variables for position and direction.
  Vector   los(1), pos(1);
  
  //--- Get complete polar angle grid
  //    1st part: the nstreams/2 Double Gauss quad angles for internal Disort use
  //    2nd part: the scat_za_grid directions (za_grid as well as
  //              cloudbox_incoming_field contains all angles of scat_za_grid,
  //              but cloudbox_incoming_field is only calculated for the
  //              downwelling ones)
  Index nn=nstreams/2;
  Index nsza = scat_za_grid.nelem();
  Index nza = nn+nsza;
  Vector za_grid(nza,0.);
  za_grid[Range(nn,nsza)] = scat_za_grid;
  cb_inc_field.resize(f_grid.nelem(),nza);
  cb_inc_field = NAN;

  Vector gmu(nn);
  Vector gwt(nn);
  
  // Call disort's gaussian quad points & weights subroutine
  qgausn_(&nn,
          gmu.get_c_array(),
          gwt.get_c_array()
         );
         
  // Calc polar angles za from their cosines mu
  for (Index i = 0; i<nn; i++)
    za_grid[i] = acos(gmu[i]) * 180./PI;

  //--- Get radiance field at boundaries (for Disort only at upper boundary)
  //    (boundary=0: lower, boundary=1: upper)
  pos[0] = z_field( cloudbox_limits[1], 0, 0 );

  for (Index za_index = 0; za_index < za_grid.nelem(); za_index ++)
    {
      los[0] =  za_grid[za_index];
      if( los[0] <=90. )
        {
          get_iy( ws, iy, t_field, z_field, vmr_field, 0, f_grid, pos, los, 
                  Vector(0), iy_unit, iy_main_agenda );

          cb_inc_field(joker,za_index) = iy(joker,0);
        }
    }
}

#else /* ENABLE_DISORT */

void get_cb_inc_field(Workspace&,
                       Matrix&,
                       const Agenda&,
                       const Tensor3&,
                       const Tensor3&,
                       const Tensor4&,
                       const Index&,
                       const ArrayOfIndex&,
                       const Vector&,
                       const Vector&,
                       const Index&
                      )
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
}

#endif /* ENABLE_DISORT */


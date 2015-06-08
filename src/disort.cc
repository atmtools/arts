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
#include "array.h"
#include "agenda_class.h"
#include "messages.h"
#include "xml_io.h"
#include "logic.h"
#include "check_input.h"
#include "auto_md.h"
#include "interpolation.h"

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
  
  \author Claudia Emde
  \date   2006-02-10
*/
void dtauc_ssalbCalc(Workspace& ws,
                     VectorView dtauc,
                     VectorView ssalb,
                     const Agenda& opt_prop_part_agenda,
                     const Agenda& propmat_clearsky_agenda,
                     const Agenda& spt_calc_agenda,
                     ConstTensor4View pnd_field,
                     ConstTensor3View t_field,
                     ConstTensor3View z_field, 
                     ConstVectorView p_grid,
                     ConstTensor4View vmr_field,
                     ConstVectorView f_mono
                    )
{
  
  const Index N_se = pnd_field.nbooks();
  // In DISORT the "cloudbox" must cover the whole atmosphere
  const Index Np_cloud = pnd_field.npages();
  const Index stokes_dim = 1; 

  assert( dtauc.nelem() == Np_cloud-1);
  assert( ssalb.nelem() == Np_cloud-1);

  // Local variables to be used in agendas
  Matrix abs_vec_spt_local(N_se, stokes_dim, 0.);
  Tensor3 ext_mat_spt_local(N_se, stokes_dim, stokes_dim, 0.);
  Matrix abs_vec_local;
  Tensor3 ext_mat_local;
  Numeric rtp_temperature_local; 
  Numeric rtp_pressure_local;
  Tensor4 propmat_clearsky_local;
  Vector ext_vector(Np_cloud); 
  Vector abs_vector(Np_cloud); 
  Vector rtp_vmr_local(vmr_field.nbooks());
  // Calculate ext_mat, abs_vec and sca_vec for all pressure points. 

  propmat_clearsky_local = 0.;
  

 for(Index scat_p_index_local = 0; scat_p_index_local < Np_cloud; 
      scat_p_index_local ++)
   {
     rtp_temperature_local = 
       t_field(scat_p_index_local, 0, 0);
     
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

 // Calculate layer averaged single scattering albedo and optical depth
 for (Index i = 0; i < Np_cloud-1; i++)
   {
     Numeric ext = 0.;
     Numeric abs = 0.;
 
     ext=.5*(ext_vector[i]+ext_vector[i+1]);
     abs=.5*(abs_vector[i]+abs_vector[i+1]);

     if (ext!=0)
       ssalb[Np_cloud-2-i]=(ext-abs)/ext;
     
     rtp_pressure_local = 0.5 * (p_grid[i] + p_grid[i+1]);
     rtp_temperature_local = 0.5 * (t_field(i,0,0) + t_field(i+1,0,0));
     
     // Average vmrs
     for (Index j = 0; j < vmr_field.nbooks(); j++)
       rtp_vmr_local[j] = 0.5 * (vmr_field(j, i, 0, 0) +
                                      vmr_field(j, i+1, 0, 0));
   
    const Vector rtp_mag_dummy(3,0);
    const Vector ppath_los_dummy;

    Tensor4 src_dummy; //FIXME: do this right
    propmat_clearsky_agendaExecute(ws,
                                  propmat_clearsky_local,
				  src_dummy,
                                  f_mono,  // monochromatic calculation
                                  rtp_mag_dummy,ppath_los_dummy,
                                  rtp_pressure_local, 
                                  rtp_temperature_local, 
                                  rtp_temperature_nlte_local_dummy,
                                  rtp_vmr_local,
                                  propmat_clearsky_agenda);  

     Numeric abs_total = propmat_clearsky_local(joker,0,0,0).sum(); //Assuming non-polarized light and only one frequency

     dtauc[Np_cloud-2-i]=(ext+abs+abs_total)*
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
  \param scat_data_mono use arts -d for docu
  \param pnd_field use arts -d for docu
  
  \author Claudia Emde
  \date   2006-02-10
*/
void phase_functionCalc(//Output
                       MatrixView phase_function,
                       //Input
                       const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                       ConstTensor4View pnd_field)
{
  Matrix phase_function_level(pnd_field.npages(), 
                              scat_data_mono[0][0].za_grid.nelem(), 0.);
  
  //Loop over pressure levels
  for (Index i_p = 0; i_p < pnd_field.npages(); i_p++)
    {
      // Loop over scattering angles
      for (Index i_t = 0; i_t < scat_data_mono[0][0].za_grid.nelem(); i_t++)
        {
          // Calculate ensemble averaged extinction coefficient
          Numeric sca_coeff=0.;

          Index i_se_flat = 0;
          for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
          {
              for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
              {
                  sca_coeff +=  pnd_field(i_se, i_p, 0, 0) *
                  (scat_data_mono[i_ss][i_se].ext_mat_data(0, 0, 0, 0, 0)-
                   scat_data_mono[i_ss][i_se].abs_vec_data(0, 0, 0, 0, 0));
                  i_se_flat++;
              }
          }

          // Phase function
          i_se_flat = 0;
          for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
          {
              for (Index i_se = 0; i_se < scat_data_mono.nelem(); i_se++)
              {
                  if (sca_coeff != 0)
                      phase_function_level(i_p, i_t) +=
                      pnd_field(i_se, i_p, 0, 0) *
                      scat_data_mono[i_ss][i_se].pha_mat_data(0, 0, i_t, 0, 0, 0, 0)
                      *4*PI/sca_coeff;// Normalization
                  i_se_flat++;
              }
          }

        }
    }


  // Calculate average phase function for the layers
  for (Index i_l = 0; i_l < pnd_field.npages()-1; i_l++)
    {
      for (Index i_t=0; i_t < phase_function_level.ncols(); i_t++)
        {
          if (phase_function_level(i_l, i_t) !=0 &&
              phase_function_level(i_l+1, i_t) !=0)
            phase_function(i_l, i_t) = .5* 
              (phase_function_level(i_l, i_t)+
               phase_function_level(i_l+1, i_t));
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
  
  \author Claudia Emde
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
        if (abs(2.-pint) > 1e-4)
        {
          CREATE_OUT1;
          out1 << "Warning: The phase function is not normalized to 2\n"
               << "The value is:" << pint << "\n";
        }
       
        pmom(i_l, joker)= 0.; 

        for (Index i = 0; i<za_grid.nelem()-1; i++) 
          {
            p0_1=1.;
            p0_2=1.;
            
            pmom(phase_function.nrows()-1-i_l,0)=1.;

            //pmom(phase_function.nrows()-1-i_l,0)+=0.5*0.5*(phase_int(i_l, i)+ 
            //                                               phase_int(i_l, i+1))
            //*abs(u[i+1]-u[i]); 
            
            p1_1=u[i];
            p1_2=u[i+1];
            
            pmom(phase_function.nrows()-1-i_l,1)+=0.5*0.5*
              (p1_1*phase_int(i_l, i)+
               p1_2*phase_int(i_l, i+1))
              *abs(u[i+1]-u[i]);
            
            for (Index l=2; l<n_legendre; l++)
              {
              p2_1=(2*(double)l-1)/(double)l*u[i]*p1_1-((double)l-1)/
                (double)l*p0_1; 
              p2_2=(2*(double)l-1)/(double)l*u[i+1]*p1_2-((double)l-1)/
                (double)l*p0_2;
              
              pmom(phase_function.nrows()-1-i_l, l)+=0.5*0.5*
                (p2_1*phase_int(i_l, i)+
                 p2_2*phase_int(i_l, i+1))
                *abs(u[i+1]-u[i]);
              
              p0_1=p1_1;
              p0_2=p1_2;
              p1_1=p2_1;
              p1_2=p2_2;
              }
          }
        // cout << "pmom : " << pmom(phase_function.nrows()-1-i_l, joker) << endl;
        
      }
    }
}

//! planck
/*! 
  Calculates the Planck function for a single temperature.
  
  Comment by CE: 
  Copied here from physics_funcs.cc, because I cannot include physics.h 
  in m_disort.cc. The problem is that the definition of complex is not 
  compatible with f2c.h
 
  Note that this expression gives the intensity for both polarisations.
  
  \return     blackbody radiation
  \param  f   frequency
  \param  t   temperature
  
  \author Patrick Eriksson 
  \date   2000-04-08 
*/
Numeric planck2( 
        const Numeric&   f, 
        const Numeric&   t )
{
  assert( f > 0 );
  assert( t >= 0 );

  // Double must be used here (if not, a becomes 0 when using float)
  static const double a = 2 * PLANCK_CONST / (SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const double b = PLANCK_CONST / BOLTZMAN_CONST;
  
  return   a * f*f*f / ( exp( b*f/t ) - 1 );
}


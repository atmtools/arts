/* Copyright (C) 2006 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                         
                           
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
  
/*!
  \file   m_disort.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   2006-02-06
  
  \brief  This file contains functions to use the multiple scattering 
  program DISORT.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "messages.h"
#include "xml_io.h"
#include "m_general.h"
#include "wsv_aux.h"
#include "disort.h"
//#include "physics_funcs.h"
#include "DISORT.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric COSMIC_BG_TEMP;


//! ScatteringDisort
/*!
  See the the online help (arts -d FUNCTION_NAME)

  \author Claudia Emde
  \date 2006-02-10
*/
void ScatteringDisort(// WS Output:
                      // WS Input
                      const ArrayOfIndex& cloudbox_limits,
                      const Index& stokes_dim,
                      const Agenda& opt_prop_part_agenda,
                      const Agenda& scalar_gas_absorption_agenda, 
                      const Agenda& spt_calc_agenda,
                      const Tensor4& pnd_field,
                      const Tensor3& t_field, 
                      const Tensor3& z_field, 
                      const Vector& p_grid, 
                      const Tensor4& vmr_field,
                      const ArrayOfSingleScatteringData& scat_data_mono,
                      const Vector& /*f_grid*/,
                      const Index& /*f_index*/, 
                      const Vector& /*scat_za_grid*/,
                      const Matrix& /*surface_emissivity_field*/)
{
  if(cloudbox_limits.nelem() != 2) 
    throw runtime_error(
                        "The cloudbox dimension is not 1D! \n" 
                        "DISORT can only be used for 1D! \n" );
  
  if (stokes_dim != 1) 
    throw runtime_error( "DISORT can only be used for unpolarized \n"
                         "calculations (i.e., stokes_dim=1),\n" );
  
  //NOTE: It is at the moment not possible to combine particle types with 
  // beinng stored on different scattering angle grids.
  // Ask whether this is required. Temperature dependance also not yet 
  // implemented. 

  const Index nlyr=cloudbox_limits[1]-cloudbox_limits[0];

  Vector dtauc(nlyr, 0.); 
  Vector ssalb(nlyr, 0.);

  dtauc_ssalbCalc(dtauc, ssalb, opt_prop_part_agenda,
                 scalar_gas_absorption_agenda, spt_calc_agenda, pnd_field, 
                 cloudbox_limits, t_field, z_field, p_grid, vmr_field);

  cout << "dtauc " << dtauc << endl
       << "ssalb " << ssalb << endl;
  
  Matrix phase_function(nlyr,scat_data_mono[0].za_grid.nelem(), 0.);
  Vector scat_angle_grid(scat_data_mono[0].za_grid.nelem(), 0.);
  scat_angle_grid = scat_data_mono[0].za_grid;
  
  phase_functionCalc(phase_function, scat_data_mono, pnd_field);
  
  cout << "phase function" << phase_function(15,joker) << "\n";  
  
  const Index n_legendre=8;
  Matrix pmom(n_legendre, nlyr, 0.); 
  pmomCalc(pmom, phase_function, scat_angle_grid, n_legendre);

  cout << "pmom " << pmom(joker, 15) << "\n";

  // Definition of other input variables for disort calculation

  /* Commented because at the moment it is not possible to convert 
     from arts Vector, Matrix to C floats, real

  const Numeric wvnmlo = f_grid[f_index]/(100*SPEED_OF_LIGHT);
  const Numeric wvnmhi = wvnmlo;
  // calculate radiant quantities at boundary of computational layers. 
  const Index usrtau = false; 
  //dummies
  const Index ntau = 0; 
  const Vector utau(1,0.);
  
  // Number of streams, I think in microwave 8 is more that sufficient
  const Index nstr=n_legendre;
  
  // Intensities to be computed for user defined polar (zenith angles)
  const Index usrang = true;
  const Index numu = scat_za_grid.nelem();
  Vector umu(numu); 
  // Transform to mu, starting with negative values
  for (Index i = 0; i<numu; i++)
    umu[i]=cos(scat_za_grid[numu-i]*PI/180);
  
  // Since we have no solar source there is no angular dependance
  const Index nphi = 1; 
  const Vector phi(nphi, 0.);
  
  const Index ibcnd=0;
  
  // Properties of solar beam, set to zero as they are not needed
  const Numeric fbeam =0.;
  const Numeric umu0=0.;
  const Numeric phi0=0.; 
  
  // Cosmic background
  const Numeric fisot = planck( f_grid[f_index], COSMIC_BG_TEMP );
  
  // surface, Lambertian if set to 0
  const Index lamber = true;
  const Numeric albedo = surface_emissivity_field(0,0);
  // only needed for bidirectional reflecting surface
  const Vector hl(1,0.); 

  //temperature of surface and cloudbox top
  const Numeric btemp = t_field(0,0,0); 
  const Numeric ttemp = t_field(cloudbox_limits[1], 0, 0); 
  
  // emissivity of top boundary, set to zero as a start, I think 
  // I have to put the gas absorption coefficient here
  const Numeric temis = 0.;
  
  // we don't need delta-scaling in microwave region
  const Index deltam = false; 
  
  // include thermal emission (very important)
  const Index plank = true; 
  
  // calculate also intensities, not only fluxes
  const Index onlyfl = false; 
  
  // Convergence criterium
  const Numeric accur = 0.005;
  
  // Specify what to be printed
  ArrayOfIndex prnt(5); 
  prnt[0]=1;
  prnt[1]=1;
  prnt[2]=1;
  prnt[3]=1;
  prnt[4]=1;
  
  const String header;
   
  const Index maxcly = 200; // Maximum number of layers
  const Index maxulv = 200; // Maximum number of user defined tau
  const Index maxumu = 200; // maximum number of zenith angles
  const Index maxcmu = 100; // maximum number of Legendre polynomials 
  const Index maxphi = 1;  //no azimuthal dependance

  // Declaration of Output variables
  Vector rfldir; 
  Vector rfldn;
  Vector flup;
  Vector dfdt;
  Vector uavg; 
  Tensor3 uu(scat_za_grid.nelem(), nlyr,1); // Intensity 
  Matrix u0u(scat_za_grid.nelem(), nlyr); // Azimuthally averaged intensity 
  Vector albmed(scat_za_grid.nelem()); // Albodo of cloudbox
  Vector trnmed(scat_za_grid.nelem()); // Transmissivity 
  
  Index header_len;

#ifdef DISORT
  disort_((integer*)nlyr, (float*) dtauc, ssalb, pmom, t_field(joker,0,0), wvnmlo, wvnmhi,
         usrtau, ntau, utau, nstr, usrang, numu, umu, nphi, phi, ibcnd, fbeam,
          umu0, phi0, fisot, lamber, albedo, hl, btemp, ttemp, temis, deltam, 
          plank, onlyfl, accur, prnt, header, maxcly, maxulv, maxumu, maxcmu, 
          maxphi, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, albmed, trnmed, 
          header_len);
#else
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
#endif
  */
}




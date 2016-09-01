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
  
/*!
  \file   m_disort.cc
  \author Claudia Emde <claudia.emde@dlr.de>
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
#include <cstring>
#include <cmath>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "messages.h"
#include "xml_io.h"
#include "m_general.h"
#include "wsv_aux.h"
//#include "disort.h"
//#include "physics_funcs.h"
//#include "disort_DISORT.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric COSMIC_BG_TEMP;



/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetFullAtm(//WS Output
                       Index& cloudbox_on,
                       ArrayOfIndex& cloudbox_limits,
                       // WS Input
                       const Vector& p_grid,
                       const Verbosity&)
{
  cloudbox_on = 1;
  cloudbox_limits.resize(2); 
  cloudbox_limits[0] = 0;
  cloudbox_limits[1] = p_grid.nelem()-1;
}
  

/*
// * Workspace method: Doxygen documentation will be auto-generated * //
#ifdef ENABLE_DISORT
void ScatteringDisort(Workspace& ws,
                      // WS Output:
                      Tensor7& scat_i_p,
                      Tensor7& scat_i_lat,
                      Tensor7& scat_i_lon,
                      Index& f_index, 
                      ArrayOfSingleScatteringData& scat_data_array_mono,
                      Tensor4& doit_i_field1D_spectrum,
                      // WS Input
                      const Index&     atmfields_checked,
                      const Index&     atmgeom_checked,
                      const Index&     cloudbox_checked,
                      const ArrayOfIndex& cloudbox_limits, 
                      const Index& stokes_dim,
                      const Agenda& opt_prop_part_agenda,
                      const Agenda& abs_scalar_gas_agenda, 
                      const Agenda& spt_calc_agenda,
                      const Tensor4& pnd_field,
                      const Tensor3& t_field, 
                      const Tensor3& z_field, 
                      const Vector& p_grid, 
                      const Tensor4& vmr_field,
                      const ArrayOfSingleScatteringData& scat_data_array,
                      const Vector& f_grid,
                      const Vector& scat_za_grid,
                      const Matrix& surface_emissivity_field,
                      const Verbosity& verbosity)
{
  CREATE_OUT1;

  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );

  out1<< "Start DISORT calculation...\n";
  
  if(pnd_field.ncols() != 1) 
    throw runtime_error(
                        "*pnd_field* is not 1D! \n" 
                        "DISORT can only be used for 1D! \n" );
  
  if (stokes_dim != 1) 
    throw runtime_error( "DISORT can only be used for unpolarized \n"
                         "calculations (i.e., stokes_dim=1),\n" );
  
  // NOTE: It is at the moment not possible to combine particle types  
  // being stored on different scattering angle grids.
  // Ask whether this is required. Temperature dependance also not yet 
  // implemented. 

  // DISORT calculations are done over the whole atmosphere because it is
  // only possible to give a constant value, i.e. cosmic background, 
  // as input, not a radiation field at the to of the domain
  Index nlyr=pnd_field.npages()-1;

  // Check whether cloudbox expands over the whole atmosphere
  
  if(cloudbox_limits.nelem()!=2 ||  cloudbox_limits[0] != 0 ||
     cloudbox_limits[1] != pnd_field.npages()-1)
    throw runtime_error("The cloudbox is not set correctly for DISORT.\n"
                        "Please use *cloudboxSetFullAtm*. \n");

  scat_i_p.resize(f_grid.nelem(), 2, 1, 1, 
                  scat_za_grid.nelem(), 1, 1);

  scat_i_p = 0.;

  doit_i_field1D_spectrum.resize(f_grid.nelem(), pnd_field.npages(), 
                                 scat_za_grid.nelem(), 1); 
 
  doit_i_field1D_spectrum= 0; 
  // Scat_i_lat, lon ---> only dummies, not used further 
  scat_i_lat.resize(1,1,1,1,1,1,1);
  scat_i_lat = 0.;
  scat_i_lon.resize(1,1,1,1,1,1,1);
  scat_i_lon = 0.; 

  // Input variables for DISORT
  // Optical depth of layers
  Vector dtauc(nlyr, 0.); 
  // Single scattering albedo of layers
  Vector ssalb(nlyr, 0.);
  
  // Phase function
  Matrix phase_function(nlyr,scat_data_array[0].za_grid.nelem(), 0.);
  // Scattering angle grid, assumed here that it is the same for
  // all particle types
  Vector scat_angle_grid(scat_data_array[0].za_grid.nelem(), 0.);
  scat_angle_grid = scat_data_array[0].za_grid;
  
  // Number of streams, I think in microwave 8 is more that sufficient
  Index nstr=8 ; 
  Index n_legendre=nstr+1;
  
  // Legendre polynomials of phase function
  Matrix pmom(nlyr, n_legendre, 0.); 

  // Intensities to be computed for user defined polar (zenith angles)
  Index usrang = TRUE_;
  Index numu = scat_za_grid.nelem();
  Vector umu(numu); 
  // Transform to mu, starting with negative values
  for (Index i = 0; i<numu; i++)
    umu[i]=cos(scat_za_grid[numu-i-1]*PI/180);
  
  // Since we have no solar source there is no angular dependance
  Index nphi = 1; 
  Vector phi(nphi, 0.);
  
  Index ibcnd=0;
  
  // Properties of solar beam, set to zero as they are not needed
  Numeric fbeam =0.;
  Numeric umu0=0.;
  Numeric phi0=0.; 

  // surface, Lambertian if set to TRUE_ 
  Index lamber = TRUE_;
  Numeric albedo = 1-surface_emissivity_field(0,0);
  // only needed for bidirectional reflecting surface
  Vector hl(1,0.); 
  
  //temperature of surface and cloudbox top
  Numeric btemp = t_field(0,0,0);
  Numeric ttemp = t_field(cloudbox_limits[1], 0, 0); 
      
  // Top of the atmosphere, emissivity = 0
  Numeric temis = 0.;
  
  // we don't need delta-scaling in microwave region
  Index deltam = FALSE_; 
      
  // include thermal emission (very important)
  Index plank = TRUE_; 
      
  // calculate also intensities, not only fluxes
  Index onlyfl = FALSE_; 
      
  // Convergence criterium
  Numeric accur = 0.005;
      
  // Specify what to be printed --> normally nothing
  Index *prnt = new Index[7]; 
  prnt[0]=FALSE_; // Input variables
  prnt[1]=FALSE_; // fluxes
  prnt[2]=FALSE_; // azimuthally averaged intensities at user 
  //and comp. angles
  prnt[3]=FALSE_; // azimuthally averaged intensities at user levels
  //and angles
  prnt[4]=FALSE_; // intensities at user levels and angles
  prnt[5]=FALSE_; // planar transmissivity and albedo 
  prnt[6]=FALSE_; // phase function moments
  
  char header[127];
  memset (header, 0, 127);
  
  Index maxcly = nlyr; // Maximum number of layers
  Index maxulv = nlyr+1; // Maximum number of user defined tau
  Index maxumu = scat_za_grid.nelem(); // maximum number of zenith angles
  Index maxcmu = n_legendre-1; // maximum number of Legendre polynomials 
  Index maxphi = 1;  //no azimuthal dependance
  
  // Declaration of Output variables
  Vector rfldir(maxulv); 
  Vector rfldn(maxulv);
  Vector flup(maxulv);
  Vector dfdt(maxulv);
  Vector uavg(maxulv);
  Tensor3 uu(maxphi, maxulv, scat_za_grid.nelem(), 0.); // Intensity 
  Matrix u0u(maxulv, scat_za_grid.nelem()); // Azimuthally averaged intensity 
  Vector albmed(scat_za_grid.nelem()); // Albedo of cloudbox
  Vector trnmed(scat_za_grid.nelem()); // Transmissivity 
      
  Vector t(nlyr+1);
 
  for (Index i = 0; i < t.nelem(); i++)
    {
      t[i] = t_field(cloudbox_limits[1]-i,0,0);
    }
      
  //dummies
  Index ntau = 0; 
  Vector utau(maxulv,0.);
  
  // Loop over frequencies
  for (f_index = 0; f_index < f_grid.nelem(); f_index ++) 
    {
      dtauc=0.;
      ssalb=0.;
      phase_function=0.;
      pmom=0.;
      
      scat_data_array_monoCalc(scat_data_array_mono, scat_data_array, f_grid, f_index, verbosity);
      
      dtauc_ssalbCalc(ws, dtauc, ssalb, opt_prop_part_agenda,
                      abs_scalar_gas_agenda, spt_calc_agenda, 
                      pnd_field, 
                      t_field, z_field, p_grid, vmr_field, f_grid[Range(f_index,1)]);
      //cout << "dtauc " << dtauc << endl
      //     << "ssalb " << ssalb << endl;
      
      phase_functionCalc(phase_function, scat_data_array_mono, pnd_field);
      //cout << "phase function" << phase_function(15,joker) << "\n";  
      
      pmomCalc(pmom, phase_function, scat_angle_grid, n_legendre, verbosity);
      //for (Index i=0; i<nlyr; i++)
      //    cout << "pmom " << pmom(i,joker) << "\n";
      
      // Wavenumber in [1/cm^2]
      Numeric wvnmlo = f_grid[f_index]/(100*SPEED_OF_LIGHT);
      Numeric wvnmhi = wvnmlo;
      
      // calculate radiant quantities at boundary of computational layers. 
      Index usrtau = FALSE_; 
      
      // Cosmic background
      Numeric fisot = planck2( f_grid[f_index], COSMIC_BG_TEMP );

        DEBUG_VAR(dtauc)
      // Call disort
      disort_(&nlyr, dtauc.get_c_array(),
              ssalb.get_c_array(), pmom.get_c_array(), 
              t.get_c_array(), &wvnmlo, &wvnmhi,
              &usrtau, &ntau, utau.get_c_array(), 
              &nstr, &usrang, &numu, 
              umu.get_c_array(), &nphi,
              phi.get_c_array(), 
              &ibcnd, &fbeam,
              &umu0, &phi0, &fisot, 
              &lamber, 
              &albedo, hl.get_c_array(),
              &btemp, &ttemp, &temis, 
              &deltam, 
              &plank, &onlyfl, &accur,
              prnt, header, 
              &maxcly, &maxulv,
              &maxumu, &maxcmu, 
              &maxphi, rfldir.get_c_array(), 
              rfldn.get_c_array(),
              flup.get_c_array(), dfdt.get_c_array(), 
              uavg.get_c_array(),
              uu.get_c_array(), u0u.get_c_array(), 
              albmed.get_c_array(),
              trnmed.get_c_array());

            cout << "intensity " << uu << endl; 
      
      
      for(Index j = 0; j<numu; j++)
        {
          for(Index k = 0; k< nlyr; k++)
            doit_i_field1D_spectrum(f_index, k+1, j, 0) = 
              uu(0,nlyr-k-1,j);
          
          scat_i_p(f_index, 0, 0, 0, j, 0, 0) = 
            uu(0, nlyr-1, j );
          scat_i_p(f_index, 1, 0, 0, j, 0, 0) = 
            uu(0, 0, j);
        }
    }
  delete [] prnt;
    
#else
void ScatteringDisort(Workspace&,
                      // WS Output:
                      Tensor7&,
                      Tensor7&,
                      Tensor7&,
                      Index&,
                      ArrayOfSingleScatteringData&,
                      Tensor4&,
                      // WS Input
                      const Index&,
                      const Index&,
                      const Index&,
                      const ArrayOfIndex&,
                      const Index&,
                      const Agenda&,
                      const Agenda&,
                      const Agenda&,
                      const Tensor4&,
                      const Tensor3&,
                      const Tensor3&,
                      const Vector&,
                      const Tensor4&,
                      const ArrayOfSingleScatteringData&,
                      const Vector&,
                      const Vector&,
                      const Matrix&,
                      const Verbosity&)
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
#endif // *ENABLE_DISORT*
}
*/ 

// No Disort support in ARTS2.2. Use ARTS2.3 or higher.
/* Workspace method: Doxygen documentation will be auto-generated */
 void ScatteringDisort(Workspace&,
                      // WS Output:
                      Tensor7&,
                      Tensor7&,
                      Tensor7&,
                      Index&,
                      ArrayOfSingleScatteringData&,
                      Tensor4&,
                      // WS Input
                      const Index&,
                      const Index&,
                      const Index&,
                      const ArrayOfIndex&,
                      const Index&,
                      const Agenda&,
                      const Agenda&,
                      const Agenda&,
                      const Tensor4&,
                      const Tensor3&,
                      const Tensor3&,
                      const Vector&,
                      const Tensor4&,
                      const ArrayOfSingleScatteringData&,
                      const Vector&,
                      const Vector&,
                      const Matrix&,
                      const Verbosity&) 
{
  throw runtime_error ("No DISORT support in ARTS2.2.\n"
                       "For using DISORT, switch to ARTS2.3 or higher.");
}


 

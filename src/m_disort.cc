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
#include "disort.h"
#include "disort_DISORT.h"
#include "messages.h"
#include "m_general.h"
#include "rte.h"
#include "wsv_aux.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric COSMIC_BG_TEMP;



/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetDisort(//WS Output
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
  

#ifdef ENABLE_DISORT
/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalc(Workspace& ws,
                // WS Output:
                Tensor7& doit_i_field,
                Index& f_index,
                ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                // WS Input
                const Index& disort_is_initialized,
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits, 
                const Agenda& propmat_clearsky_agenda, 
                const Agenda& opt_prop_part_agenda,
                const Agenda& spt_calc_agenda,
                const Agenda& iy_main_agenda,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Vector& scat_za_grid,
                const Vector& surface_scalar_reflectivity,
                const Index& nstreams,
                const Index& non_iso_inc,
                const Verbosity& verbosity)
{
  CREATE_OUT1;
  CREATE_OUT0;

  // NOTE: It is at the moment not possible to combine scattering elements 
  // being stored on different scattering angle grids. Ask if this is required.
  // Temperature dependence also not yet implemented. 

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) return;

  // Check whether DisortInit was executed
  if (!disort_is_initialized)
    {
      ostringstream os;
      os << "Initialization method *DisortInit* must be called before "
         << "*DisortCalc*";
      throw runtime_error( os.str() );
    }

  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );

  if( pnd_field.ncols() != 1 ) 
    throw runtime_error(
                        "*pnd_field* is not 1D! \n" 
                        "DISORT can only be used for 1D! \n" );
  
  // DISORT requires even number of streams:
  // nstreams is total number of directions, up- and downwelling, and the up-
  // and downwelling directions need to be symmetrically distributed, i.e. same
  // number of directions in both hemispheres is required. horizontal direction
  // (90deg) can not be covered in a plane-parallel atmosphere.
  if( nstreams/2*2 != nstreams )
    {
      ostringstream os;
      os << "DISORT requires an even number of streams, but yours is "
         << nstreams << ".\n";
      throw runtime_error( os.str() );
    }

  if( surface_scalar_reflectivity.nelem() != f_grid.nelem()  &&  
      surface_scalar_reflectivity.nelem() != 1 )
    {
      ostringstream os;
      os << "The number of elements in *surface_scalar_reflectivity* should\n"
         << "match length of *f_grid* or be 1."
         << "\n length of *f_grid* : " << f_grid.nelem() 
         << "\n length of *surface_scalar_reflectivity* : " 
         << surface_scalar_reflectivity.nelem()
         << "\n";
      throw runtime_error( os.str() );
    }

  if( min(surface_scalar_reflectivity) < 0  ||  
      max(surface_scalar_reflectivity) > 1 )
    {
      throw runtime_error( 
         "All values in *surface_scalar_reflectivity* must be inside [0,1]." );
    }

  doit_i_field.resize(f_grid.nelem(), pnd_field.npages(), 1, 1,
                      scat_za_grid.nelem(), 1, 1);
  doit_i_field = 0.;

  // Input variables for DISORT
  Index nlyr;
  bool do_non_iso;
  if( p_grid.nelem() == pnd_field.npages() )
    // cloudbox covers whole atmo anyways. no need for a calculation of
    // non-iso incoming field at top-of-cloudbox. disort will be run over
    // whole atmo.
    {
      do_non_iso = false;
      nlyr = p_grid.nelem()-1;
    }
  else
    {
      if( non_iso_inc )
        // cloudbox only covers part of atmo. disort will be initialized with
        // non-isotropic incoming field and run over cloudbox only.
        {
          do_non_iso = true;
          nlyr = pnd_field.npages()-1;
        }
      else
        // cloudbox only covers part of atmo. disort will be run over whole
        // atmo, though (but only in-cloudbox rad field passed to doit_i_field).
        {
          do_non_iso = false;
          nlyr = p_grid.nelem()-1;
        }
    }
      
  // Optical depth of layers
  Vector dtauc(nlyr, 0.); 
  // Single scattering albedo of layers
  Vector ssalb(nlyr, 0.);
  
  // Phase function
  Matrix phase_function(nlyr,scat_data[0][0].za_grid.nelem(), 0.);
  // Scattering angle grid, assumed here that it is the same for
  // all scattering elements
  Vector scat_angle_grid(scat_data[0][0].za_grid.nelem(), 0.);
  scat_angle_grid = scat_data[0][0].za_grid;
  
  Index nstr=nstreams;
  Index n_legendre=nstreams+1;
  
  // Legendre polynomials of phase function
  Matrix pmom(nlyr, n_legendre, 0.); 

  // Intensities to be computed for user defined polar (zenith angles)
  Index usrang = TRUE_;
  Index numu = scat_za_grid.nelem();
  Vector umu(numu); 
  // Transform to mu, starting with negative values
  for (Index i = 0; i<numu; i++)
    umu[i] = -cos(scat_za_grid[i]*PI/180);

  
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
  // only needed for bidirectional reflecting surface
  Vector hl(1,0.); 
  // albedo only set in freq-loop (as it might be freq-dependent
  
  //temperature of surface
  Numeric btemp = t_field(0,0,0);

  //upper boundary conditions:
  // DISORT offers isotropic incoming radiance or emissivity-scaled planck
  // emission. Both are applied additively.
  // We want to have cosmic background radiation, for which ttemp=COSMIC_BG_TEMP
  // and temis=1 should give identical results to fisot(COSMIC_BG_TEMP). As they
  // are additive we should use either the one or the other.
  // Note: previous setup (using fisot) setting temis=0 should be avoided.
  // Generally, temis!=1 should be avoided since that technically implies a
  // reflective upper boundary (though it seems that this is not exploited in
  // DISORT1.2, which we so far use).

  // Cosmic background
  // we use temis*ttemp as upper boundary specification, hence CBR set to 0.
  Numeric fisot = 0;

  // Top of the atmosphere temperature and emissivity
  //Numeric ttemp = t_field(cloudbox_limits[1], 0, 0); 
  //Numeric temis = 0.;
  Numeric ttemp = COSMIC_BG_TEMP;
  Numeric temis = 1.;

  // Top of the atmosphere non-isotropic incoming radiation
  Matrix cb_inc_field;
  Vector intang(scat_za_grid.nelem()+nstr/2, 0.);
  if( do_non_iso )
    get_cb_inc_field( ws, cb_inc_field,
                      iy_main_agenda,
                      z_field, t_field, vmr_field, cloudbox_limits,
                      f_grid, scat_za_grid, nstreams );

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
  if( nstr<4 )
    maxcmu = 4; // reset for low nstr since DISORT selftest uses 4 streams,
                // hence requires at least 4 Legendre polynomials
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
      t[i] = t_field(nlyr-i,0,0);
  
  //dummies
  Index ntau = 0; 
  Vector utau(maxulv,0.);
  
  // Loop over frequencies
  for (f_index = 0; f_index < f_grid.nelem(); f_index ++) 
    {
      Numeric albedo;
      if( surface_scalar_reflectivity.nelem()>1 )
        albedo = surface_scalar_reflectivity[f_index];
      else
        albedo = surface_scalar_reflectivity[0];
      
      // Top of the atmosphere non-isotropic incoming radiation
      if( do_non_iso )
        {
          // extract monchromatic field from cloudbox_incoming_field
          intang = cb_inc_field(f_index,joker);

          // Moved this assert into DISORT.f. There we can test the actually
          // applied intang values for validity (and skip upwelling angle values
          // at the end of intang, which are deliberately set to NaN.).
          //for( Index i_za=0; i_za<intang.nelem(); i_za++ )
          //  assert( !(isnan(intang[i_za]) || intang[i_za]<0.) );

          // convert ARTS units to DISORT units
          // W/(m2 sr Hz) -> W/(m2 sr cm-1)
          intang *= (100*SPEED_OF_LIGHT);
          // we replace the isotropic TOA source by the non-isotropic incoming
          // one, hence set TOA source (via source temperature) to 0
          ttemp = 0.;
        }
      else
        {
          intang = 0.;
          ttemp = COSMIC_BG_TEMP;
        }

      scat_data_monoCalc(scat_data_mono, scat_data, f_grid, f_index, verbosity);
      
      dtauc_ssalbCalc(ws, dtauc, ssalb,
                      propmat_clearsky_agenda,
                      spt_calc_agenda, opt_prop_part_agenda,
                      pnd_field,
                      t_field(Range(0,nlyr+1),joker,joker),
                      z_field(Range(0,nlyr+1),joker,joker),
                      vmr_field(joker,Range(0,nlyr+1),joker,joker),
                      p_grid[Range(0,nlyr+1)],
                      cloudbox_limits, f_grid[Range(f_index,1)]);

      phase_functionCalc(phase_function, scat_data_mono, pnd_field,
                         cloudbox_limits);

      for( Index l=0; l<nlyr; l++ )
        if( phase_function(l,0)==0. )
          assert( ssalb[l]==0. );

      pmomCalc(pmom, phase_function, scat_angle_grid, n_legendre, verbosity);
      
      // Wavenumber in [1/cm]
      Numeric wvnmlo = f_grid[f_index]/(100*SPEED_OF_LIGHT);
      Numeric wvnmhi = wvnmlo;
      
      // calculate radiant quantities at boundary of computational layers. 
      Index usrtau = FALSE_; 
      
      //DEBUG_VAR(dtauc)
      
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
              intang.get_c_array(),
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

      for(Index j = 0; j<numu; j++)
          for(Index k = 0; k<(cloudbox_limits[1]-cloudbox_limits[0]+1); k++)
            doit_i_field(f_index, k, 0, 0, j, 0, 0) =
              uu(0,nlyr-k-cloudbox_limits[0],j) / (100*SPEED_OF_LIGHT);
    }
  delete [] prnt;

}

#else /* ENABLE_DISORT */

void DisortCalc(Workspace&,
                // WS Output:
                Tensor7&,
                Index&,
                ArrayOfArrayOfSingleScatteringData&,
                // WS Input
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const ArrayOfIndex&,
                const Agenda&,
                const Agenda&,
                const Agenda&,
                const Agenda&,
                const Tensor4&,
                const Tensor3&,
                const Tensor3&,
                const Tensor4&,
                const Vector&,
                const ArrayOfArrayOfSingleScatteringData&,
                const Vector&,
                const Vector&,
                const Vector&,
                const Index&,
                const Index&,
                const Verbosity&)
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
}

#endif /* ENABLE_DISORT */


/* Workspace method: Doxygen documentation will be auto-generated */
void DisortInit(//WS Output
              Tensor7& doit_i_field,
              Index& disort_is_initialized,
              // WS Input
              const Index& stokes_dim,
              const Index& atmosphere_dim,
              const Vector& f_grid,
              const Vector& scat_za_grid,
              const Index& cloudbox_on,
              const ArrayOfIndex& cloudbox_limits,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Verbosity& verbosity)
{
  if (!cloudbox_on)
  {
    CREATE_OUT0;
    disort_is_initialized = 0;
    out0 << "  Cloudbox is off, scattering calculation will be skipped.\n";
    return;
  }
  
  // -------------- Check the input ------------------------------
  
  if( atmosphere_dim != 1   )
    throw runtime_error( "For running DISORT, atmospheric dimensionality "
                         "must be 1.\n");

  if (stokes_dim < 0 || stokes_dim > 1)
    throw runtime_error( "For running DISORT, the dimension of stokes vector "
                         "must be 1.\n");

  if( cloudbox_limits[0] != 0   )
    {
      ostringstream os;
      os << "DISORT calculations currently only possible with\n"
         << "lower cloudbox at 0th atmospheric level.\n";
      throw runtime_error(os.str());
    }

  // Zenith angle grid.
  Index nza = scat_za_grid.nelem();

  if (scat_za_grid[0] != 0. || scat_za_grid[nza-1] != 180.)
    throw runtime_error( "The range of *scat_za_grid* must [0 180]." );
  
  if (!is_increasing(scat_za_grid))
    throw runtime_error("*scat_za_grid* must be increasing.");

  // scat_za_grid here is only relevant to provide an i_field from which the
  // sensor los angles can be interpolated by yCalc; it does not the determine
  // the accuracy of the DISORT output itself at these angles. So we can only
  // apply a very rough test here, whether the grid is appropriate. However, we
  // set the threshold fairly high since calculation costs for a higher number
  // of angles are negligible.
  if ( nza < 37 )
    {
      ostringstream os;
      os << "We require size of scat_za_grid to be > 36\n"
         << "to ensure accurate radiance field interpolation in yCalc.\n"
         << "Note that for DISORT additional computation costs for\n"
         << "larger numbers of angles are negligible.";
      throw runtime_error( os.str() );
    }

  if( nza/2*2 != nza )
    {
      // uneven nza detected. uneven nza (when set as equidistant grid as
      // commonly done by ARTS) lead to polar angle grid point at 90deg, ie at
      // the horizontal. this is not safely calculable in a plane-parallel atmo.
      // for now we just force the user to use an even nza.
      //
      // an even nza does not place the center angles close to horizon, though,
      // unless the number of streams is very high. therefore, one could instead
      // replace this gridpoint with two points centered closely around 90deg
      // and derive the 90deg value from averaging these two.
      // however, this is left to the future (and needs testing).
      //
      // FIXME: more correct (and stable in case of non-equidistant grids) is to
      // check whether scat_za_grid actually contains the 90deg angle and to
      // reject (or circumvent) this specifically.
      ostringstream os;
      os << "Uneven nza detected. nza=" << nza << ".\n";
      throw runtime_error( os.str() );
    }

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  if ( scat_data.empty() )
    throw runtime_error(
                         "No single scattering data present.\n"
                         "See documentation of WSV *scat_data* for options to "
                         "make single scattering data available.\n"
                         );

  // DISORT can only handle randomly oriented particles.
  bool allp20=true;
  for( Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++ )
    for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
      if( scat_data[i_ss][i_se].ptype != PTYPE_MACROS_ISO )
        allp20=false;
  if( !allp20 )
    {
      ostringstream os;
      os << "DISORT can only handle scattering elements of type "
         << PTYPE_MACROS_ISO << " (" << PTypeToString(PTYPE_MACROS_ISO) << "),\n"
         << "but at least one element of other type (" << PTYPE_HORIZ_AL
         << "=" << PTypeToString(PTYPE_HORIZ_AL) << " or " << PTYPE_GENERAL
         << "=" << PTypeToString(PTYPE_GENERAL) << ") is present.\n";
      throw runtime_error( os.str() );
    }
    
  // So far our interface can only handle particles with single scattering data
  // given on identical angular grids.
  const Vector data_za_grid = scat_data[0][0].za_grid;
  const Index ndza = data_za_grid.nelem();
  bool ident_anggrid=true;
  for( Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++ )
    for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
      // not an exhaustive test, but should catch most cases: checking identical
      // size as well as identical second and second to last elements. no use in
      // checking first and last elements as they should be 0 and 180 and this
      // should have been checked elsewhere.
      if( scat_data[i_ss][i_se].za_grid.nelem() != ndza ||
          scat_data[i_ss][i_se].za_grid[1] != data_za_grid[1] ||
          scat_data[i_ss][i_se].za_grid[ndza-2]!=data_za_grid[ndza-2] )
        ident_anggrid=false;
   if( !ident_anggrid )
    {
      ostringstream os;
      os << "ARTS-DISORT currently requires identical angular grids of\n"
         << "scattering data for all scattering elements, but yours differ.\n";
      throw runtime_error( os.str() );
    }

  //------------- end of checks ---------------------------------------
  
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index Nza = scat_za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  doit_i_field.resize( Nf, Np_cloud, 1, 1, Nza, 1, 1 );
  doit_i_field = NAN;
  
  disort_is_initialized = 1;
}


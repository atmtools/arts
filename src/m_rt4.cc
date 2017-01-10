/* Copyright (C) 2016 Jana Mendrok <jana.mendrok@gmail.com>

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
  \file   m_rt4.cc
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Workspace functions related to application of scattering solver RT4.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "auto_md.h"
#include <complex.h>
#include "messages.h"
#include "m_xml.h"
#include "rt4.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric COSMIC_BG_TEMP;

#ifdef ENABLE_RT4
/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Calc( Workspace& ws,
                // WS Output:
                Tensor7& doit_i_field,
                Vector& scat_za_grid,
                Vector& scat_aa_grid,
                Index& f_index,
                ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                // WS Input
                const Index& rt4_is_initialized,
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits, 
                const Agenda& propmat_clearsky_agenda, 
                const Agenda& opt_prop_part_agenda,
                const Agenda& spt_calc_agenda,
                //const Agenda& iy_main_agenda,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Index& stokes_dim,
                //const Vector& scat_za_grid,
                const Vector& surface_scalar_reflectivity,
                const Index& nstreams,
                const Index& non_iso_inc,
                const String& pfct_method,
                const String& quad_type,
                const Index& pfct_aa_grid_size,
                const Numeric& pfct_threshold,
                const Numeric& max_delta_tau,
                const Verbosity& verbosity )
{
  CREATE_OUT1;
  CREATE_OUT0;

  // NOTE: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) return;

  // Check whether DisortInit was executed
  if (!rt4_is_initialized)
    {
      ostringstream os;
      os << "Initialization method *RT4Init* must be called before "
         << "*RT4Calc*";
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
    throw runtime_error("*pnd_field* is not 1D! \n" 
                        "RT4 can only be used for 1D!\n" );

  Index nhza;
  if( quad_type=="D" || quad_type=="G" )
    {
      nhza=1;
    }
  else if( quad_type=="L" )
    {
      nhza=0;
    }
  else
    {
      ostringstream os;
      os << "Unknown quadrature type.\n";
      throw runtime_error(os.str());
    }

  //cout << "doit_i_field has " << doit_i_field.npages() << "angles.\n";
  //cout << "considering " << nstreams << " streams and " << 2*nhza
  //     << " extra angles.\n";
  if( doit_i_field.npages() != nstreams+2*nhza ) 
    throw runtime_error("Sizes of *doit_i_field* is inconsistent with *nstreams*.\n"
                        "Make sure to use the same *nstreams* (and *quad_type*)"
                        " with *RT4Init* and *RT4Calc!\n" );
  
  // RT4 actually uses number of angles in single hemisphere. However, we don't
  // want a bunch of different approaches used in the interface, so we apply the
  // DISORT (and ARTS) way of total number of angles here. Hence, we have to
  // ensure here that total number is even.
  if( nstreams/2*2 != nstreams )
    {
      ostringstream os;
      os << "RT4 requires an even number of streams, but yours is "
         << nstreams << ".\n";
      throw runtime_error( os.str() );
    }
  Index nhstreams=nstreams/2;
  // nummu is the total number of angles in one hemisphere, including the
  // quadrature angles as well as the "extra" angles.
  Index nummu=nhstreams+nhza;

  if( pfct_method!="interpolate" )
  {
    // The old interface can only handle particles with single scattering data
    // given on identical angular grids.
    const Vector data_za_grid = scat_data[0][0].za_grid;
    const Index ndza = data_za_grid.nelem();
    bool ident_anggrid=true;
    for( Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++ )
      for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
        // not an exhaustive test, but should catch most cases: checking
        // identical size as well as identical second and second to last
        // elements. no use in checking first and last elements as they should
        // be 0 and 180 and this should have been checked elsewhere.
        if( scat_data[i_ss][i_se].za_grid.nelem() != ndza ||
            scat_data[i_ss][i_se].za_grid[1] != data_za_grid[1] ||
            scat_data[i_ss][i_se].za_grid[ndza-2]!=data_za_grid[ndza-2] )
          ident_anggrid=false;
     if( !ident_anggrid )
      {
        ostringstream os;
        os << "ARTS-RT4 currently requires identical angular grids of\n"
           << "scattering data for all scattering elements, but yours differ.\n";
        throw runtime_error( os.str() );
      }
  }

  // FIXME: remove/replace the following two tests when other than
  // ground_type=="L" is implemented.
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

  // Input variables for RT4
  Index num_layers=p_grid.nelem()-1;
  /*
  bool do_non_iso;
  if( p_grid.nelem() == pnd_field.npages() )
    // cloudbox covers whole atmo anyways. no need for a calculation of
    // non-iso incoming field at top-of-cloudbox. disort will be run over
    // whole atmo.
    {
      do_non_iso = false;
      num_layers = p_grid.nelem()-1;
    }
  else
    {
      if( non_iso_inc )
        // cloudbox only covers part of atmo. disort will be initialized with
        // non-isotropic incoming field and run over cloudbox only.
        {
          do_non_iso = true;
          num_layers = pnd_field.npages()-1;
        }
      else
        // cloudbox only covers part of atmo. disort will be run over whole
        // atmo, though (but only in-cloudbox rad field passed to doit_i_field).
        {
          do_non_iso = false;
          num_layers = p_grid.nelem()-1;
        }
    }
  */
      
  const Index nstokes=stokes_dim;

  // temperature of surface
  const Numeric ground_temp = t_field(0,0,0);

  // surface reflection type.
  //  proprietary RT4 allows (L)ambertian and (F)resnel.
  //  FIXME: For first, I hardcode that to Lambertian. Implementation of Fresnel
  //  as well as specular reflection with given reflectivity to be done later.
  //  For that review handling of theseoptions in ARTS (yCalc) and do that
  //  consistently here (e.g. applying agendas...).
  const String ground_type="L";

  // Top of the atmosphere temperature
  //  FIXME: so far hard-coded to cosmic background. However, change that to set
  //  according to/from space_agenda.
  const Numeric sky_temp = COSMIC_BG_TEMP;

  // Data fields
  Vector height(num_layers+1);
  Vector temperatures(num_layers+1);
  for (Index i = 0; i < height.nelem(); i++)
    {
      height[i] = z_field(num_layers-i,0,0);
      temperatures[i] = t_field(num_layers-i,0,0);
    }

  // this indexes all cloudbox layers as cloudy layers.
  // optional FIXME: to use the power of RT4 (faster solving scheme for
  // individual non-cloudy layers), one should consider non-cloudy layers within
  // cloudbox. That requires some kind of recognition and index setting based on
  // pnd_field or (derived) cloud layer extinction or scattering.
  const Index num_scatlayers=pnd_field.npages()-1;
  Vector scatlayers(num_layers, 0.);
  Vector gas_extinct(num_layers, 0.);
  Tensor6 scatter_matrix(num_scatlayers,4,nummu,nstokes,nummu,nstokes, 0.);
  Tensor5 extinct_matrix(num_scatlayers,2,nummu,nstokes,nstokes, 0.);
  Tensor4 emis_vector(num_scatlayers,2,nummu,nstokes, 0.);

  for (Index i = 0; i < cloudbox_limits[1]-cloudbox_limits[0]; i++)
    {
      scatlayers[num_layers-1-cloudbox_limits[0]-i] = float(i+1);
    }

  // in RT4 mu_values is generally only output. however, we need the values for
  // preparing the single scattering data at these angles. therefore, we
  // calculate them here using RT4's proprietary quadrature methods. They
  // simultaneously provide the quadrature weights, too. We keep them so far,
  // might use them for ensuring proper normalization in the preparation of the
  // single scattering data.
  Vector mu_values(nummu, 0.);
  Vector quad_weights(nummu, 0.);

  if( quad_type=="D" )
    {
      double_gauss_quadrature_(nhstreams,
                               mu_values.get_c_array(),
                               quad_weights.get_c_array()
                              );
    }
  else if( quad_type=="G" )
    {
      gauss_legendre_quadrature_(nhstreams,
                                 mu_values.get_c_array(),
                                 quad_weights.get_c_array()
                                );
    }
  else //if( quad_type=="L" )
    {
      lobatto_quadrature_(nhstreams,
                          mu_values.get_c_array(),
                          quad_weights.get_c_array()
                         );
    }

  // Set "extra" angle (at 0deg) if quad_type!="L"
  if( nhza>0 )
      mu_values[nhstreams] = 1.;

  // FIXME: we should be able to avoid setting scat_za_grid here in one way,
  // and resetting in another before leaving the WSM. This, however, requires
  // rearranging the angle order and angle assignment in the RT4-SSP prep
  // routines.
  scat_za_grid.resize(2*nummu);
  for (Index imu=0; imu<nummu; imu++)
    {
      scat_za_grid[imu] = acos(mu_values[imu]) * RAD2DEG;
      scat_za_grid[nummu+imu] = 180.-scat_za_grid[imu];
      //cout << "Setting za[" << imu << "]=" << scat_za_grid[imu]
      //     << " and  za[" << nummu+imu << "]=" << scat_za_grid[nummu+imu]
      //     << " from mu[" << imu << "]=" << mu_values[imu] << "\n";
    }
  scat_aa_grid.resize(1);
  scat_aa_grid[0] = 0.;

  // Output variables
  Tensor3 up_rad(num_layers+1,nummu,nstokes, 0.);
  Tensor3 down_rad(num_layers+1,nummu,nstokes, 0.);


  // Loop over frequencies
  for (f_index = 0; f_index < f_grid.nelem(); f_index ++) 
    {
      // surface albedo
      //  (used only if surface assumed to be Lambertian)
      Numeric ground_albedo;
      if( surface_scalar_reflectivity.nelem()>1 )
        ground_albedo = surface_scalar_reflectivity[f_index];
      else
        ground_albedo = surface_scalar_reflectivity[0];

      // surface refractive index
      //  (only used in case of Fresnel surface. As Fresnel not yet implemented,
      //  this is hardcoded as a dummy so far).
      Complex ground_index;
      ground_index.real(-1.);
      ground_index.imag(-1.);

      // Wavelength [um]
      Numeric wavelength;
      wavelength = 1e6*SPEED_OF_LIGHT/f_grid[f_index];
      //cout << "# processing freq #" << f_index << " at " << f_grid[f_index]*1e-9
      //     << "GHz\n";

      scat_data_monoCalc(scat_data_mono, scat_data, f_grid, f_index, verbosity);
      
      gas_optpropCalc( ws, gas_extinct,
                       propmat_clearsky_agenda,
                       t_field(Range(0,num_layers+1),joker,joker),
                       vmr_field(joker,Range(0,num_layers+1),joker,joker),
                       p_grid[Range(0,num_layers+1)],
                       f_grid[Range(f_index,1)]);

      par_optpropCalc( ws, emis_vector, extinct_matrix,
                       //scatlayers,
                       spt_calc_agenda, opt_prop_part_agenda,
                       pnd_field,
                       t_field(Range(0,num_layers+1),joker,joker),
                       cloudbox_limits, stokes_dim, nummu );
      sca_optpropCalc( scatter_matrix,
                       emis_vector, extinct_matrix,
                       scat_data_mono, pnd_field, stokes_dim,
                       scat_za_grid, quad_weights,
                       pfct_method, pfct_aa_grid_size, pfct_threshold,
                       verbosity );

//#pragma omp critical(fortran_rt4)
//      {
          // Call RT4
          radtrano_(nstokes,
               nummu,
               nhza,
               max_delta_tau,
               quad_type.c_str(),
               ground_temp,
               ground_type.c_str(),
               ground_albedo,
               ground_index,
               sky_temp,
               wavelength,
               num_layers,
               height.get_c_array(),
               temperatures.get_c_array(),
               gas_extinct.get_c_array(),
               num_scatlayers,
               scatlayers.get_c_array(),
               extinct_matrix.get_c_array(),
               emis_vector.get_c_array(),
               scatter_matrix.get_c_array(),
               //noutlevels,
               //outlevels.get_c_array(),
               mu_values.get_c_array(),
               up_rad.get_c_array(),
               down_rad.get_c_array()
                 );
//      }

      // RT4 rad output is in wavelength units, nominally in W/(m2 sr um), where
      // wavelength input is required in um.
      // FIXME: When using wavelength input in m, output should be in W/(m2 sr
      // m). However, check this. So, at first we use wavelength in um. Then
      // change and compare.
      //    
      // FIXME: if ever we allow the cloudbox to be not directly at the surface
      // (at atm level #0, respectively), the assigning from up/down_rad to
      // doit_i_field needs to checked. there seems some offsetting going on
      // (test example: TestDOIT.arts. if kept like below, doit_i_field at
      // top-of-cloudbox seems to actually be from somewhere within the
      // cloud(box) indicated by downwelling being to high and downwelling
      // exhibiting a non-zero polarisation signature (which it wouldn't with
      // only scalar gas abs above).
      //
      Numeric rad_l2f = wavelength/f_grid[f_index];
      // down/up_rad contain the radiances in order from slant (90deg) to steep
      // (0 and 180deg, respectively) streams,then the possible extra angle(s).
      // We need to resort them properly into doit_i_field, such that order is
      // from 0 to 180deg.
      for(Index j = 0; j<nummu; j++)
        for(Index k = 0; k<(cloudbox_limits[1]-cloudbox_limits[0]+1); k++)
          for (Index ist = 0; ist<stokes_dim; ist++ )
              {
                //doit_i_field(f_index, k, 0, 0, nummu+j, 0, ist) =
                //  up_rad(num_layers-k,j,ist)*rad_l2f;
                //doit_i_field(f_index, k, 0, 0, nummu-1-j, 0, ist) =
                //  down_rad(num_layers-k,j,ist)*rad_l2f;
                doit_i_field(f_index, k, 0, 0, nummu+j, 0, ist) =
                  up_rad(num_layers-k,j,ist)*rad_l2f;
                doit_i_field(f_index, k, 0, 0, nummu-1-j, 0, ist) =
                  down_rad(num_layers-k,j,ist)*rad_l2f;
              }
    }

  //scat_za_grid.resize(nstreams);
  //cout << "Setting scat_za_grid for ARTS consistent with doit_i_field.\n";
  for (Index j=0; j<nummu; j++)
    {
      scat_za_grid[nummu-1-j] = acos(mu_values[j])*RAD2DEG;
      scat_za_grid[nummu+j] = 180.-acos(mu_values[j])*RAD2DEG;
      //cout << "setting scat_za[" << nummu-1-j << "]=" << scat_za_grid[nummu-1-j]
      //     << " and [" << nummu+j << "]=" << scat_za_grid[nummu+j]
      //     << " from mu[" << j << "]=" << mu_values[j] << "\n";
    }
}

#else /* ENABLE_RT4 */

void RT4Calc( Workspace&,
                // WS Output:
                Tensor7&,
                Vector&,
                Vector&,
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
                const Tensor4&,
                const Tensor3&, 
                const Tensor3&, 
                const Tensor4&,
                const Vector&, 
                const ArrayOfArrayOfSingleScatteringData&,
                const Vector&,
                const Index&,
                //const Vector&,
                const Vector&,
                const Index&,
                const Index&,
                const String&,
                const String&,
                const Index&,
                const Numeric&,
                const Numeric&,
                const Verbosity& )
{
    throw runtime_error ("This version of ARTS was compiled without RT4 support.");
}

#endif /* ENABLE_RT4 */


/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Init(//WS Output
              Tensor7& doit_i_field,
              Index& rt4_is_initialized,
              // WS Input
              const Index& stokes_dim,
              const Index& atmosphere_dim,
              const Vector& f_grid,
              const Index& cloudbox_on,
              const ArrayOfIndex& cloudbox_limits,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Index& nstreams,
              const String& quad_type,
              const Verbosity& verbosity)
{
  if (!cloudbox_on)
  {
    CREATE_OUT0;
    rt4_is_initialized = 0;
    out0 << "  Cloudbox is off, scattering calculation will be skipped.\n";
    return;
  }
  
  // -------------- Check the input ------------------------------
  
  if( atmosphere_dim != 1   )
    throw runtime_error( "For running RT4, atmospheric dimensionality "
                         "must be 1.\n");

  if (stokes_dim < 0 || stokes_dim > 2)
    throw runtime_error( "For running RT4, the dimension of stokes vector "
                         "must be 1 or 2.\n");

  if( cloudbox_limits[0] != 0   )
    {
      ostringstream os;
      os << "RT4 calculations currently only possible with "
         << "lower cloudbox limit\n"
         << "at 0th atmospheric level "
         << "(assumes surface there, ignoring z_surface).\n";
      throw runtime_error(os.str());
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

  if( quad_type.length()>1 )
    {
      ostringstream os;
      os << "Input parameter *quad_type* not allowed to contain more than a "
         << "single character.\n"
         << "Yours has " << quad_type.length() << ".\n";
      throw runtime_error(os.str());
    }

  // If quad!='L' we need to add extra angle 0 & 180deg for final scat_za_grid
  Index neza;
  if( quad_type=="D" || quad_type=="G" )
    {
      neza=2;
    }
  else if( quad_type=="L" )
    {
      neza=0;
    }
  else
    {
      ostringstream os;
      os << "Unknown quadrature type.\n";
      throw runtime_error(os.str());
    }

  // RT4 can only completely or azimuthally randomly oriented particles.
  bool no_p10=true;
  for( Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++ )
    for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
      if( scat_data[i_ss][i_se].ptype != PTYPE_MACROS_ISO &&
          scat_data[i_ss][i_se].ptype != PTYPE_HORIZ_AL )
        no_p10=false;
  if( !no_p10 )
    {
      ostringstream os;
      os << "RT4 can only handle scattering elements of type "
         << PTYPE_MACROS_ISO << " (" << PTypeToString(PTYPE_MACROS_ISO) << ") and\n"
         << PTYPE_HORIZ_AL << " (" << PTypeToString(PTYPE_HORIZ_AL) << "),\n"
         << "but at least one element of other type (" << PTYPE_GENERAL
         << "=" << PTypeToString(PTYPE_GENERAL) << ") is present.\n";
      throw runtime_error( os.str() );
    }
    
  //------------- end of checks ---------------------------------------
  
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  //const Index Nza = scat_za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  //doit_i_field.resize( Nf, Np_cloud, 1, 1, Nza, 1, 1 );
  doit_i_field.resize( Nf, Np_cloud, 1, 1, nstreams+neza, 1, stokes_dim );
  doit_i_field = NAN;
  
  rt4_is_initialized = 1;
}



/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef ENABLE_RT4
void RT4Test( Tensor4& out_rad,
              const String& datapath,
              const Verbosity& verbosity )
{
    rt4_test( out_rad, datapath, verbosity );
}
#else
void RT4Test( Tensor4&,
              const String&,
              const Verbosity& )
{
    throw runtime_error ("This version of ARTS was compiled without RT4 support.");
}
#endif


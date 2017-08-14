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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "agenda_class.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "disort.h"
#include "disort_DISORT.h"
#include "interpolation.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "rte.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric COSMIC_BG_TEMP;

//! check_disort_input
/*!
  Checks that input of DisortCalc* is sane.

  \param cloudbox_on           as the WSV 
  \param disort_is_initialized  as the WSV 
  \param atmfields_checked     as the WSV 
  \param atmgeom_checked       as the WSV 
  \param cloudbox_checked      as the WSV 
  \param scat_data             as the WSV
  \param scat_za_grid          as the WSV
  \param nstreams              Number of quadrature angles (both hemispheres).
  \param pfct_method           see DisortCalc doc.
  \param pnd_ncols             Number of columns (latitude points) in *pnd_field*.
  \param ifield_npages         Number of pages (polar angle points) in *doit_i_field*.
  
  \author Jana Mendrok
  \date   2017-02-23
*/
void check_disort_input( // Input
                         const Index& cloudbox_on,
                         const Index& atmfields_checked,
                         const Index& atmgeom_checked,
                         const Index& cloudbox_checked,
                         const Index& scat_data_checked,
                         const Index& atmosphere_dim,
                         const Index& stokes_dim,
                         const ArrayOfIndex& cloudbox_limits, 
                         const ArrayOfArrayOfSingleScatteringData& scat_data,
                         ConstVectorView scat_za_grid,
                         const Index& nstreams,
                         const String& pfct_method,
                         const Index& pnd_ncols )
{
  // Don't do anything if there's no cloudbox defined.
  //if (!cloudbox_on) return;
  // Seems to loopholy to just skip the scattering, so rather throw an error
  // (assuming if RT4 is called than it's expected that a scattering calc is
  // performed. semi-quietly skipping can easily be missed and lead to wrong
  // conclusions.).
  if (!cloudbox_on)
  {
    throw runtime_error( "Cloudbox is off, no scattering calculations to be"
                         "performed." );
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
  if( scat_data_checked != 1 )
    throw runtime_error( "The scat_data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );

  if( atmosphere_dim != 1   )
    throw runtime_error( "For running DISORT, atmospheric dimensionality "
                         "must be 1.\n");

  if (stokes_dim < 0 || stokes_dim > 1)
    throw runtime_error( "For running DISORT, the dimension of stokes vector "
                         "must be 1.\n");

  if( pnd_ncols != 1 ) 
    throw runtime_error("*pnd_field* is not 1D! \n" 
                        "DISORT can only be used for 1D!\n" );

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  if( cloudbox_limits[0] != 0   )
    {
      ostringstream os;
      os << "DISORT calculations currently only possible with "
         << "lower cloudbox limit\n"
         << "at 0th atmospheric level "
         << "(assumes surface there, ignoring z_surface).\n";
      throw runtime_error(os.str());
    }

  if ( scat_data.empty() )
    throw runtime_error(
                         "No single scattering data present.\n"
                         "See documentation of WSV *scat_data* for options to "
                         "make single scattering data available.\n"
                         );

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

  // Zenith angle grid.
  Index nza = scat_za_grid.nelem();

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

  if (scat_za_grid[0] != 0. || scat_za_grid[nza-1] != 180.)
    throw runtime_error( "The range of *scat_za_grid* must [0 180]." );
  
  if (!is_increasing(scat_za_grid))
    throw runtime_error("*scat_za_grid* must be increasing.");

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
      os << "Uneven number of angles in scat_za_grid (nza=" << nza << ".\n"
         << "Use even number with no grid point at 90deg poin.t\n";
      throw runtime_error( os.str() );
    }

  // DISORT can only handle randomly oriented particles.
  bool all_p20=true;
  for( Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++ )
    for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
      if( scat_data[i_ss][i_se].ptype != PTYPE_TOTAL_RND )
        all_p20=false;
  if( !all_p20 )
    {
      ostringstream os;
      os << "DISORT can only handle scattering elements of type "
         << PTYPE_TOTAL_RND << " (" << PTypeToString(PTYPE_TOTAL_RND) << "),\n"
         << "but at least one element of other type (" << PTYPE_AZIMUTH_RND
         << "=" << PTypeToString(PTYPE_AZIMUTH_RND) << " or " << PTYPE_GENERAL
         << "=" << PTypeToString(PTYPE_GENERAL) << ") is present.\n";
      throw runtime_error( os.str() );
    }
    
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
        os << "ARTS-DISORT currently requires identical angular grids of\n"
           << "scattering data for all scattering elements, but yours differ.\n";
        throw runtime_error( os.str() );
      }
  }
}


//! init_ifield
/*!
  Initialize doit_i_field with the right size and NaN values.

  \param doit_i_field          as the WSV
  \param f_grid                as the WSV
  \param cloudbox_limits       as the WSV
  \param nang                  Total number of angles with RT output.
  \param stokes_dim            as the WSV

  \author Jana Mendrok
  \date   2017-03-06
*/
void init_ifield( // Output
                  Tensor7& doit_i_field,
                  // Input
                  const Vector& f_grid,
                  const ArrayOfIndex& cloudbox_limits, 
                  const Index& nang,
                  const Index& stokes_dim )
{
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  //const Index Nza = scat_za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  //doit_i_field.resize( Nf, Np_cloud, 1, 1, Nza, 1, 1 );
  doit_i_field.resize( Nf, Np_cloud, 1, 1, nang, 1, stokes_dim );
  doit_i_field = NAN;
}


//! get_disortsurf_props
/*!
  Derive surface property input for RT4's proprietary surface handling depending
  on surface reflection type.

  \param ground_albedo         Scalar surface albedo (for ground_type=L).
  \param ground_reflec         Vector surface relfectivity (for ground_type=S).
  \param ground_index          Surface complex refractive index (for ground_type=F).
  \param f_grid                as the WSV
  \param ground_type           Surface reflection type flag.
  \param surface_skin_t        as the WSV
  \param surface_scalar_reflectivity  as the WSV (used with ground_type=L).
  \param surface_reflectivity  as the WSV (used with ground_type=S).
  \param surface_complex_refr_index   as the WSV (used with ground_type=F).
  \param stokes_dim            as the WSV

  \author Jana Mendrok
  \date   2017-02-23
*/
void get_disortsurf_props( // Output
                        Vector& albedo,
                        Numeric& btemp,
                        // Input
                        ConstVectorView f_grid,
                        const Numeric& surface_skin_t,
                        ConstVectorView surface_scalar_reflectivity )
{
  // temperature of surface
  if (surface_skin_t<0. || surface_skin_t>1000.)
    {
      ostringstream os;
      os << "Surface temperature has been set or derived as " << btemp << " K,\n"
         << "which is not considered a meaningful value.\n"
         << "For surface method 'L', *surface_skin_t* needs to\n"
         << "be set and passed explicitly. Maybe you didn't do this?";
      throw runtime_error( os.str() );
    }
  btemp = surface_skin_t;

  // surface albodo
  if( surface_scalar_reflectivity.nelem() != f_grid.nelem()  &&  
      surface_scalar_reflectivity.nelem() != 1 )
    {
      ostringstream os;
      os << "The number of elements in *surface_scalar_reflectivity*\n"
         << "should match length of *f_grid* or be 1."
         << "\n length of *f_grid* : " << f_grid.nelem() 
         << "\n length of *surface_scalar_reflectivity* : " 
         << surface_scalar_reflectivity.nelem()
         << "\n";
      throw runtime_error( os.str() );
    }

  if ( min(surface_scalar_reflectivity) < 0  ||  
       max(surface_scalar_reflectivity) > 1 )
    {
      throw runtime_error( "All values in *surface_scalar_reflectivity*"
                           " must be inside [0,1]." );
    }

  if( surface_scalar_reflectivity.nelem()>1 )
    for (Index f_index = 0; f_index < f_grid.nelem(); f_index ++) 
      albedo[f_index] = surface_scalar_reflectivity[f_index];
  else
    for (Index f_index = 0; f_index < f_grid.nelem(); f_index ++) 
      albedo[f_index] = surface_scalar_reflectivity[0];
}

//! dtauc_ssalbCalc
/*!
  Calculates layer averaged cloud optical depth (dtauc) and 
  single scattering albedo (ssalb). These variables are required as
  input for the DISORT subroutine

  \param ws                    Current workspace
  \param dtauc                 optical depths for all layers
  \param ssalb                 single scattering albedos for all layers
  \param scat_data             as the WSV
  \param f_index               index of frequency grid point handeled
  \param propmat_clearsky_agenda as the WSA
  \param pnd_field             as the WSV 
  \param t_field               as the WSV 
  \param z_field               as the WSV 
  \param vmr_field             as the WSV 
  \param p_grid                as the WSV 
  \param cloudbox_limits       as the WSV 
  \param f_mono                frequency (single entry vector)
  
  \author Claudia Emde, Jana Mendrok
  \date   2006-02-10
*/
void dtauc_ssalbCalc( Workspace &ws,
                      VectorView dtauc,
                      VectorView ssalb,
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      const Index& f_index,
                      const Agenda& propmat_clearsky_agenda,
                      ConstTensor4View pnd_field,
                      ConstTensor3View t_field,
                      ConstTensor3View z_field, 
                      ConstTensor4View vmr_field,
                      ConstVectorView p_grid,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstVectorView f_mono,
                      const Verbosity& verbosity )
{
  // Initialization
  dtauc=0.;
  ssalb=0.;
  
  const Index N_se = pnd_field.nbooks();

  const Index Np_cloud = pnd_field.npages();
  const Index Np = p_grid.nelem();

  assert( dtauc.nelem() == Np-1);
  assert( ssalb.nelem() == Np-1);

  const Index stokes_dim = 1; 

  // Local variables to be used in agendas
  ArrayOfStokesVector abs_vec_spt_local(N_se);
  for(auto& av : abs_vec_spt_local)
  {
    av = StokesVector(1, stokes_dim);
    av.SetZero();
  }
  
  ArrayOfPropagationMatrix ext_mat_spt_local(N_se);
  for(auto& pm : ext_mat_spt_local)
  {
    pm = PropagationMatrix(1, stokes_dim);
    pm.SetZero();
  }
    
  StokesVector abs_vec_local;
  PropagationMatrix ext_mat_local;
  Numeric rtp_temperature_local; 
  Numeric rtp_pressure_local;
  ArrayOfPropagationMatrix propmat_clearsky_local;
  Vector ext_vector(Np, 0.); 
  Vector abs_vector(Np, 0.); 
  Vector rtp_vmr_local(vmr_field.nbooks());
  Vector za_dummy(1, 0.), aa_dummy(1, 0.);

  // Calculate ext_mat and abs_vec for all pressure points in cloudbox 
  for(Index scat_p_index_local = 0;
            scat_p_index_local < Np_cloud; 
            scat_p_index_local ++)
    {
      rtp_temperature_local =
        t_field(scat_p_index_local+cloudbox_limits[0], 0, 0);
     
      //Calculate optical properties for all individual scattering elements:
      opt_prop_sptFromScat_data(ext_mat_spt_local, abs_vec_spt_local,
                                scat_data, 1,
                                za_dummy, aa_dummy, 0, 0, f_index,
                                rtp_temperature_local, pnd_field, 
                                scat_p_index_local, 0, 0, verbosity);

      opt_prop_bulkCalc(ext_mat_local, abs_vec_local, 
                        ext_mat_spt_local, abs_vec_spt_local,
                        pnd_field,
                        scat_p_index_local, 0, 0, verbosity);

      ext_vector[scat_p_index_local+cloudbox_limits[0]] = ext_mat_local.Kjj()[0];
      abs_vector[scat_p_index_local+cloudbox_limits[0]] = abs_vec_local.Kjj()[0];
    }

  const Vector  rtp_temperature_nlte_local_dummy(0);

  // Calculate layer averaged single scattering albedo and layer optical depth
  for(auto& pm : propmat_clearsky_local)
  {
    pm.SetZero();
  }
  
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
      ArrayOfStokesVector nlte_dummy;
      // This is right since there should be only clearsky partials
      ArrayOfPropagationMatrix partial_dummy;
      ArrayOfStokesVector partial_source_dummy,partial_nlte_dummy;
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
      Numeric abs_gas = 0.;
      for(auto& pm : propmat_clearsky_local)
        abs_gas += pm.Kjj()[0];

      if (ext_part!=0)
        ssalb[Np-2-i]=(ext_part-abs_part) / (ext_part+abs_gas);
     
      dtauc[Np-2-i]=(ext_part+abs_gas)*
        (z_field(i+1, 0, 0)-z_field(i, 0, 0));
    }  
}

//! phase_functionCalc2
/*!
  Calculates layer averaged normalized phase functions from 
  the phase matrix in SingleScatteringData.
  Temperature and angle grid interpolations are applied.

  \param phase_function  normalized layer-averaged bulk phase function
  \param scat_data         as the WSV
  \param f_index           index of frequency grid point handeled
  \param pnd_field         as the WSV
  \param t_field           as the WSV 
  \param cloudbox_limits   as the WSV
  \param pfct_za_grid_size number of equidistant scatt. angles in 0-180deg

  \author Claudia Emde, Jana Mendrok
  \date   2006-02-10
*/
void phase_functionCalc2( //Output
                          MatrixView phase_function,
                          //Input
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const Index& f_index,
                          ConstTensor4View pnd_field,
                          ConstTensor3View t_field,
                          const ArrayOfIndex& cloudbox_limits,
                          const Index& pfct_za_grid_size,
                          const Verbosity& verbosity )
{
  // Initialization
  phase_function=0.;

  const Index nlyr = phase_function.nrows();
  const Index Np_cloud = pnd_field.npages();
  const Index N_se = pnd_field.nbooks();
  const Index stokes_dim = 1; 

  Matrix phase_function_level(Np_cloud, 
                              pfct_za_grid_size, 0.);
  Vector sca_coeff_level(Np_cloud, 0.);

  // Local variables to be used in agendas
  Numeric rtp_temperature_local; 
  ArrayOfStokesVector abs_vec_spt_local(N_se);
  for(auto& av : abs_vec_spt_local)
  {
    av = StokesVector(1, stokes_dim);
    av.SetZero();
  }
  
  StokesVector abs_vec_local;
  ArrayOfPropagationMatrix ext_mat_spt_local(N_se);
  for(auto& pm : ext_mat_spt_local)
    pm = PropagationMatrix(1, stokes_dim);
  
  PropagationMatrix ext_mat_local;
  Tensor5 pha_mat_spt_local(N_se, pfct_za_grid_size, 1,
                            stokes_dim, stokes_dim, 0.);
  Tensor4 pha_mat_local;

  Vector za_grid;
  nlinspace(za_grid, 0, 180, pfct_za_grid_size);
  Vector aa_grid(1, 0.);

  for(Index scat_p_index_local = 0;
            scat_p_index_local < Np_cloud; 
            scat_p_index_local ++)
    {
/*
      if( scat_p_index_local < 1 )
      {
      cout << "at cloud-lev #" << scat_p_index_local << "\n";
      Index i_se_flat=0;
      Numeric sca_coeff;
      for( Index i_ss=0; i_ss<scat_data.nelem(); i_ss++ )
      {
        cout << " scat species #" << i_ss << "\n";
        for( Index i_se=0; i_se<scat_data[i_ss].nelem(); i_se++ )
        {
          if( i_se_flat==65 )
          {
          cout << "  scat element #" << i_se << " (flat element #" << i_se_flat
               << ")\n";
          for( Index i_t=0; i_t<scat_data[i_ss][i_se].T_grid.nelem(); i_t++ )
          {
            sca_coeff = scat_data[i_ss][i_se].ext_mat_data(0,i_t,0,0,0) -
                        scat_data[i_ss][i_se].abs_vec_data(0,i_t,0,0,0);
            cout << "   T[" << i_t << "]=" << scat_data[i_ss][i_se].T_grid[i_t]
                 << "K with sca_coeff=" << sca_coeff << "\n";
            if( sca_coeff!=0. )
            {
              Numeric intP=0.;
              for (Index i_a = 0; i_a <
                   scat_data[i_ss][i_se].za_grid.nelem(); i_a++)
              {
                //cout << "    processing ang grid point #" << i_a << "\n";
                if( i_a>0 )
                {
                  intP += PI *
                        (scat_data[i_ss][i_se].pha_mat_data(0,i_t,i_a,0,0,0,0)+
                         scat_data[i_ss][i_se].pha_mat_data(0,i_t,i_a,0,0,0,0)) *
                        (abs(cos(scat_data[i_ss][i_se].za_grid[i_a]*PI/180.)-
                             cos(scat_data[i_ss][i_se].za_grid[i_a-1]*PI/180.)));
                }
              }
              cout << "  at T[" << i_t << "]=" << scat_data[i_ss][i_se].T_grid[i_t]
                   << "K: int PFCT / scatcoef="
                     << intP / sca_coeff << "\n";
            }
          }
          }
          i_se_flat++;
        }
      }
      }
*/

      // Calculate scat_coef from ext_mat and abs_vec for all pressure points in
      // cloudbox 
      // FIXME: This is a copy of what is done in dtauc_ssalbCalc. Might be more
      // clever to merge these two methods (or to output scat_coef profile from
      // dtauc_ssalbCalc & input here) to avoid redundant calculations.
      rtp_temperature_local =
        t_field(scat_p_index_local+cloudbox_limits[0], 0, 0);
     
      //Calculate optical properties for all individual scattering elements:
      opt_prop_sptFromScat_data(ext_mat_spt_local, abs_vec_spt_local,
                                scat_data, 1,
                                za_grid, aa_grid, 0, 0, f_index,
                                rtp_temperature_local, pnd_field, 
                                scat_p_index_local, 0, 0, verbosity);

      opt_prop_bulkCalc(ext_mat_local, abs_vec_local, 
                        ext_mat_spt_local, abs_vec_spt_local,
                        pnd_field,
                        scat_p_index_local, 0, 0, verbosity);

      //cout << "  extmat_total=" << ext_mat_local(0,0,0) << "\n";
      //cout << "  absvec_total=" << abs_vec_local(0,0) << "\n";
      sca_coeff_level[scat_p_index_local] =
        ext_mat_local(0,0,0) - abs_vec_local(0,0);
      //cout << "  => scatcoef_total=" << sca_coeff_level[scat_p_index_local] << "\n";

      if (sca_coeff_level[scat_p_index_local] != 0)
        {
          // Calculate the phase matrix of individual scattering elements
          pha_mat_sptFromScat_data(pha_mat_spt_local,
                                   scat_data, 1,
                                   za_grid, aa_grid, 0, 0, // angles, only needed for za=0
                                   f_index, rtp_temperature_local,  pnd_field,
                                   scat_p_index_local, 0, 0, verbosity );
              
          // Sum over all scattering elements
          pha_matCalc( pha_mat_local,
                       pha_mat_spt_local, pnd_field, 1,
                       scat_p_index_local, 0, 0,
                       verbosity );

          // Bulk scattering function
          // (conversion to phase function only done when doing layer averaging.
          // this because averaging needs to be on scat coeff weighted phase
          // function aka bulk scattering function)
          phase_function_level(scat_p_index_local, joker) =
            pha_mat_local(joker, 0, 0, 0);

/*
          Numeric intP=0.;
          Vector intP_se(N_se,0.);

          for (Index i_t = 0; i_t < phase_function_level.ncols(); i_t++)
            {
              if( i_t>0 )
              {
                  intP += PI *
                    (phase_function_level(scat_p_index_local, i_t) +
                     phase_function_level(scat_p_index_local, i_t-1)) *
                    abs(cos(za_grid[i_t]*PI/180.)-
                        cos(za_grid[i_t-1]*PI/180.));
                  for( Index i_se=0; i_se<N_se; i_se++ )
                    intP_se[i_se] += PI *
                      (pha_mat_spt_local(i_se,i_t,0,0,0) +
                       pha_mat_spt_local(i_se,i_t-1,0,0,0) ) *
                      abs(cos(za_grid[i_t]*PI/180.)-
                         cos(za_grid[i_t-1]*PI/180.));
              }
            }

          if( scat_p_index_local < 1 )
          {
            cout << "at lev_cloud #" << scat_p_index_local << " (T="
                 << rtp_temperature_local << "K): "
                 << "  integ PFCT=" << intP << ", total scatcoef="
                 << sca_coeff_level[scat_p_index_local]
                 << " => integ PFCT/total scatcoef="
                 << intP/sca_coeff_level[scat_p_index_local] << "\n";
            Index i_se=65;
            cout << "  at scatelem #" << i_se << ": int PFCT="
                 << intP_se[i_se] << ", scatcoeff="
                 << ext_mat_spt_local(i_se,0,0)-abs_vec_spt_local(i_se,0)
                 << " => int PFCT/scatcoef="
                 << intP_se[i_se] / (ext_mat_spt_local(i_se,0,0)-abs_vec_spt_local(i_se,0))
                 << "\n";
*/
/*
              for( Index i_se=0; i_se<N_se; i_se++ )
              if( ext_mat_spt_local(i_se,0,0)-abs_vec_spt_local(i_se,0)!=0. &&
                intP_se[i_se]!=0. )
                {
                  cout << "  at scatelem #" << i_se << ": int PFCT / scatcoef="
                       << intP_se[i_se] /
                          (ext_mat_spt_local(i_se,0,0)-abs_vec_spt_local(i_se,0))
                       << "\n";
                }
              else
                {
                  cout << "  at scatelem #" << i_se << ": int PFCT="
                       << intP_se[i_se] << ", scatcoeff="
                       << ext_mat_spt_local(i_se,0,0)-abs_vec_spt_local(i_se,0)
                       << "\n";
                }

            }
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
          //Numeric intP=0.;
          for (Index i_t=0; i_t < phase_function_level.ncols(); i_t++)
          {
              phase_function(lyr_id, i_t) = 4*PI *
                ( phase_function_level(i_l, i_t) + 
                  phase_function_level(i_l+1, i_t) ) /
                ( sca_coeff_level[i_l] + sca_coeff_level[i_l+1] );
              //if( i_t>0 )
              //    intP += 0.5 *
              //      (phase_function(lyr_id, i_t) + phase_function(lyr_id, i_t-1)) *
              //      abs(cos(za_grid[i_t]*PI/180.)-cos(za_grid[i_t-1]*PI/180.));
          }
          //cout << "at lyr #" << lyr_id << " (from cloud levs #" << i_l
          //     << " and " << i_l+1 << "): intP=" << intP << "\n";
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


//! phase_functionCalc
/*!
  Calculates layer averaged normalized phase functions from 
  the phase matrix in SingleScatteringData. The scattering angle 
  grid is taken from the data. 
  It is required that all scattering elements are given on the same 
  scattering angle grid. No temperature interpolation done.

  \param phase_function normalized phase function
  \param scat_data       as the WSV
  \param f_index         index of frequency grid point handeled
  \param pnd_field       as the WSV
  \param cloudbox_limits as the WSV
  
  \author Claudia Emde, Jana Mendrok
  \date   2006-02-10
*/
void phase_functionCalc(//Output
                        MatrixView phase_function,
                        //Input
                        const ArrayOfArrayOfSingleScatteringData& scat_data,
                        const Index& f_index,
                        ConstTensor4View pnd_field,
                        const ArrayOfIndex& cloudbox_limits,
                        const String pfct_method)
{
  // Turned out we get some numerical issues in converting scattering matrix to
  // phase function (sca.coef vs. 4Pi normalized) if pnd's are too low. Hence,
  // set a threshold below which pnd is assumed as zero.
  // 1e-99 means less than one particle in our galaxy.
  Numeric pnd_threshold = 1e-99;

  // Initialization
  phase_function=0.;
  const Index nlyr = phase_function.nrows();

  const Index Np_cloud = pnd_field.npages();
  Matrix phase_function_level(Np_cloud, scat_data[0][0].za_grid.nelem(), 0.);
  
  Vector sca_coeff_level(Np_cloud, 0.);
  Index this_f_index;

  //Loop over pressure levels
  for (Index i_p = 0; i_p < Np_cloud; i_p++)
    {
      // Calculate ensemble averaged scattering coefficient
      Numeric sca_coeff=0.;
      Index i_se_flat=0;
      //Numeric intP=0.;

      for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
        {
          for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++)
            {
              if( pnd_field(i_se_flat, i_p, 0, 0) > pnd_threshold )
                {
                  // FIXME: In case we allow K,a,Z to have different
                  // f-dimensions, we need to derive this_f_index individually
                  // for K, a, and Z!
                  if( scat_data[i_ss][i_se].ext_mat_data.nshelves()==1 )
                    this_f_index = 0;
                  else
                    this_f_index = f_index;
    
                  // For T, we assume that K and a have the same T dimensions
                  // (but Z can have a different one). That is as checked by
                  // scat_data_checkedCalc.
                  Index this_T_index = -1;
                  if ( scat_data[i_ss][i_se].ext_mat_data.nbooks()==1 )
                    {
                      this_T_index = 0;
                    }
                  else
                    {
                      if( pfct_method=="low" )
                        this_T_index = 0;
                      else if( pfct_method=="high" )
                        this_T_index =
                          scat_data[i_ss][i_se].ext_mat_data.nbooks()-1;
                      else //if( pfct_method=="median" )
                        this_T_index =
                          scat_data[i_ss][i_se].ext_mat_data.nbooks()/2;
                    }

                  sca_coeff +=  pnd_field(i_se_flat, i_p, 0, 0) *
                    ( scat_data[i_ss][i_se].ext_mat_data(this_f_index,
                                                         this_T_index,
                                                         0, 0, 0)
                     -scat_data[i_ss][i_se].abs_vec_data(this_f_index,
                                                         this_T_index,
                                                         0, 0, 0));
                }
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
          for (Index i_t = 0; i_t < scat_data[0][0].za_grid.nelem(); i_t++)
            {
              i_se_flat=0;
              for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
              {
                for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++)
                {
                  if( scat_data[i_ss][i_se].pha_mat_data.nlibraries()==1 )
                    this_f_index = 0;
                  else
                    this_f_index = f_index;
    
                  Index this_T_index = -1;
                  if ( scat_data[i_ss][i_se].pha_mat_data.nvitrines()==1 )
                    {
                      this_T_index = 0;
                    }
                  else
                    {
                      if( pfct_method=="low" )
                        this_T_index = 0;
                      else if( pfct_method=="high" )
                        this_T_index =
                          scat_data[i_ss][i_se].pha_mat_data.nvitrines()-1;
                      else //if( pfct_method=="median" )
                        this_T_index =
                          scat_data[i_ss][i_se].pha_mat_data.nvitrines()/2;
                    }

                  phase_function_level(i_p, i_t) += 
                    pnd_field(i_se_flat, i_p, 0, 0) *
                    scat_data[i_ss][i_se].pha_mat_data(this_f_index,
                                                       this_T_index,
                                                       i_t, 0, 0, 0, 0);
                  i_se_flat++;
                }
              }
            }
/*
              if( i_t>0 )
                  intP += PI *
                    (phase_function_level(i_p, i_t) +
                     phase_function_level(i_p, i_t-1)) *
                    abs(cos(scat_data[0][0].za_grid[i_t]*PI/180.)-
                        cos(scat_data[0][0].za_grid[i_t-1]*PI/180.));
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

//! pmomCalc2
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
void pmomCalc2(//Output
              MatrixView pmom,
              //Input
              ConstMatrixView phase_function, 
              ConstVectorView scat_angle_grid,
              const Index n_legendre,
              const Verbosity& verbosity)
{
  assert( phase_function.ncols() == scat_angle_grid.nelem() );

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
             << "more than 10%. Happens for cloudbox layer #" << i_l << ".\n"
             << "Something is wrong with your scattering data. Check!\n";
          throw runtime_error( os.str() );
        }
        if (abs(2.-pint) > 2e-2)
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
  assert( phase_function.ncols() == scat_angle_grid.nelem() );

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
      //cout << "at layer #" << i_l << " P_int=" << pint << "\n";
      
      if (pint != 0){
        if (abs(2.-pint) > 0.2)
        {
          ostringstream os;
          os << "Phase function norm in layer " << i_l << " is " << pint
             << ", i.e. deviates\n"
             << "from expected value (2.0) by " << abs(2.-pint)*50. << "%.\n"
             << "Something is wrong with your scattering data. Check!\n";
          throw runtime_error( os.str() );
        }
        if (abs(2.-pint) > 2e-2)
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


//! surf_albedoCalc
/*!
  Calculates the diffuse power reflection coefficient (an estimate of the total
  surface albedo, equivalent to ARTS' surface_scalar_reflectivity) from
  reflection matrices according to *surface_rt_prop_agenda* settings for use as
  input parameter albedo to a Disort calculation (internally applying a
  Lambertian surface).

  \param albedo                Diffuse power reflection coefficient.
  \param surface_rtprop_agenda as the WSA
  \param f_grid                as the WSV
  \param scat_za_grid          as the WSV
  \param surf_alt              Surface altitude.
  
  \author Jana Mendrok
  \date   2017-02-16
*/
void surf_albedoCalc( Workspace& ws, 
                      //Output
                      VectorView albedo,
                      Numeric& btemp,
                      //Input
                      const Agenda& surface_rtprop_agenda,
                      ConstVectorView f_grid,
                      ConstVectorView scat_za_grid,
                      const Numeric& surf_alt,
                      const Verbosity& verbosity )
{
  // Here, we derive an average surface albedo of the setup as given by ARTS'
  // surface_rtprop_agenda to use with Disorts's proprietary Lambertian surface.
  // In this way, ARTS-Disort can approximately mimick all surface reflection
  // types that ARTS itself can handle.
  // Surface temperature as derived from surface_rtprop_agenda is also returned.
  //
  // We derive the reflection matrices over all incident and reflected polar
  // angle directions and derive their integrated value (or weighted average).
  // The surface_rtprop_agenda handles one reflected direction at a time
  // (rtp_los) and for the reflected directions we loop over all (upwelling)
  // angles as given by scat_za_grid. The derived reflection matrices already
  // represent the reflectivity (or power reflection coefficient) for the given
  // reflected direction including proper angle weighting. For integrating/
  // averaging over the reflected directions, we use the same approach, i.e.
  // weight each angle by its associated range as given by the half distances to
  // the neighboring grid angles (using ARTS' scat_za_grid means 0/180deg are
  // grid points, 90deg shouldn't be among them (resulting from even number
  // requirement for Disort and the (implicitly assumed?) requirement of a
  // horizon-symmetric za_grid)).
  // 
  // We do all frequencies here at once (assuming this is the faster variant as
  // the agenda anyway (can) provide output for full f_grid at once and as we
  // have to apply the same inter/extrapolation to all the frequencies).
  //
  // FIXME: Allow surface to be elsewhere than at lowest atm level (this
  // requires changes in the surface setting part and more extensive ones in the
  // atmospheric optical property prep part within the frequency loop further
  // below).

  CREATE_OUT0;

  chk_not_empty( "surface_rtprop_agenda", surface_rtprop_agenda );

  const Index nf = f_grid.nelem();
  Index frza=0;
  while( frza<scat_za_grid.nelem() && scat_za_grid[frza]<90.)
      frza++;
  if( frza==scat_za_grid.nelem() )
    {
      ostringstream os;
      os << "No upwelling direction found in scat_za_grid.\n";
      throw runtime_error(os.str());
    }
  const Index nrza=scat_za_grid.nelem()-frza;
  //cout << nrza << " upwelling directions found, starting from element #"
  //     << frza << " of scat_za_grid.\n";
  Matrix dir_refl_coeff(nrza,nf, 0.);

  // Local input of surface_rtprop_agenda.
  Vector rtp_pos(1, surf_alt); //atmosphere_dim is 1

  // first derive the (reflected-)direction dependent power reflection
  // coefficient
  for (Index rza=0; rza<nrza; rza++)
    {
      // Local output of surface_rtprop_agenda.
      Numeric   surface_skin_t;
      Matrix    surface_los;
      Tensor4   surface_rmatrix;
      Matrix    surface_emission;

      Vector rtp_los(1, scat_za_grid[rza+frza]);
      out0 << "Doing reflected dir #" << rza << " at " << rtp_los[0] << " degs\n";

      surface_rtprop_agendaExecute( ws,
                                    surface_skin_t, surface_emission,
                                    surface_los, surface_rmatrix,
                                    f_grid, rtp_pos, rtp_los,
                                    surface_rtprop_agenda );
      //cout << "surf_los has " << surface_los.ncols() << " columns and "
      //     << surface_los.nrows() << " rows.\n";
      assert( surface_los.ncols()==1 || surface_los.nrows()==0 );
      if( rza==0 )
        btemp = surface_skin_t;
      else if( surface_skin_t != btemp )
        {
          ostringstream os;
          os << "Something went wrong.\n"
             << "  *surface_rtprop_agenda* returned different surface_skin_t\n"
             << "  for different LOS.\n";
          throw runtime_error(os.str());
        }
      if(  surface_los.nrows()>0 )
        {
          for (Index f_index=0; f_index<nf; f_index++)
            dir_refl_coeff(rza,f_index) =
              surface_rmatrix(joker,f_index,0,0).sum();
        }
      out0 << "  directional albedos[f_grid] = " <<  dir_refl_coeff(rza,joker)
           << "\n";
    }

  if (btemp<0. || btemp>1000.)
    {
      ostringstream os;
      os << "Surface temperature has been derived as " << btemp << " K,\n"
         << "which is not considered a meaningful value.\n";
      throw runtime_error( os.str() );
    }

  // now integrate/average the (reflected-)direction dependent power reflection
  // coefficients
  //
  // starting with deriving the angles defining the angle ranges
  Vector surf_int_grid(nrza+1);
  // the first angle grid point should be around (but above) 90deg and should
  // cover the angle range between the 90deg and half-way point towards the next
  // angle grid point. we probably also want to check, that we don't
  // 'extrapolate' too much.
  if( is_same_within_epsilon(scat_za_grid[frza],90.,1e-6) )
    {
      ostringstream os;
      os << "Looks like scat_za_grid contains the 90deg direction,\n"
         << "which it shouldn't for running Disort.\n";
      throw runtime_error(os.str());
    }
  Numeric za_extrapol = (scat_za_grid[frza]-90.) /
                        (scat_za_grid[frza+1]-scat_za_grid[frza]);
  const Numeric ok_extrapol=0.5;
  if( (za_extrapol-ok_extrapol)>1e-6 )
    {
      ostringstream os;
      os << "Extrapolation range from shallowest scat_za_grid point\n"
         << "to horizon is too big.\n"
         << "  Allowed extrapolation factor is 0.5.\n  Yours is "
         << za_extrapol << ", which is " << za_extrapol-0.5 << " too big.\n";
      throw runtime_error(os.str());
    }
  if( !is_same_within_epsilon(scat_za_grid[scat_za_grid.nelem()-1],180.,1e-6) )
    {
      ostringstream os;
      os << "Looks like last point in scat_za_grid is not 180deg.\n";
      throw runtime_error(os.str());
    }

  surf_int_grid[0] = 90.;
  surf_int_grid[nrza] = 180.;
  for (Index rza=1; rza<nrza; rza++)
    surf_int_grid[rza] = 0.5*(scat_za_grid[frza+rza-1]+scat_za_grid[frza+rza]);
  surf_int_grid *= DEG2RAD;

  // now calculating the actual weights and apply them
  for (Index rza=0; rza<nrza; rza++)
    {
      //Numeric coslow = cos(2.*surf_int_grid[rza]);
      //Numeric coshigh = cos(2.*surf_int_grid[rza+1]);
      //Numeric w = 0.5*(coshigh-coslow);
      Numeric w = 0.5*(cos(2.*surf_int_grid[rza+1])-cos(2.*surf_int_grid[rza]));
      //cout << "at reflLOS[" << rza << "]=" << scat_za_grid[frza+rza] << ":\n";
      //cout << "  angle weight derives as w = 0.5*(" << coshigh << "-"
      //     << coslow << ") = " << w << "\n";
      //cout << "  weighting directional reflection coefficient from "
      //     <<  dir_refl_coeff(rza,joker);
      dir_refl_coeff(rza,joker) *= w;
      //cout << " to " <<  dir_refl_coeff(rza,joker) << "\n";
    }

  // eventually sum up the weighted directional power reflection coefficients
  for (Index f_index=0; f_index<nf; f_index++)
    {
      albedo[f_index] = dir_refl_coeff(joker,f_index).sum();
      out0 << "at f=" << f_grid[f_index]*1e-9 << " GHz, ending up with albedo="
           << albedo[f_index] << "\n";
      if( albedo[f_index]<0 || albedo[f_index]>1. )
        {
          ostringstream os;
          os << "Something went wrong: Albedo must be inside [0,1],\n"
             << "  but is not at freq #" << f_index << " , where it is "
             << albedo[f_index] << ".\n";
          throw runtime_error( os.str() );
        }
    }
}


#ifdef ENABLE_DISORT
//! run_disort
/*!
  Prepares actual input variables for Disort, runs it, and sorts the output into
  doit_i_field.

  \param ws                    Current workspace
  \param doit_i_field          as the WSV
  \param f_grid                as the WSV
  \param p_grid                as the WSV
  \param z_field               as the WSV
  \param t_field               as the WSV
  \param vmr_field             as the WSV
  \param pnd_field             as the WSV
  \param scat_data             as the WSV
  \param propmat_clearsky_agenda  as the WSA
  \param iy_main_agenda        as the WSA
  \param cloudbox_limits       as the WSV 
  \param surface_skin_t        as the WSV
  \param surface_scalar_reflectivity  as the WSM
  \param scat_za_grid          as the WSV
  \param nstreams              Number of quadrature angles (both hemispheres).
  \param non_iso_inc           see DisortCalc doc.
  \param pfct_method           see DisortCalc doc.

  \author Jana Mendrok
  \date   2017-02-23
*/
void run_disort( Workspace& ws,
              // Output
              Tensor7& doit_i_field,
              // Input
              ConstVectorView f_grid,
              ConstVectorView p_grid,
              ConstTensor3View z_field,
              ConstTensor3View t_field,
              ConstTensor4View vmr_field,
              ConstTensor4View pnd_field,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Agenda& propmat_clearsky_agenda, 
              const Agenda& iy_main_agenda,
              const ArrayOfIndex& cloudbox_limits,
              Numeric& surface_skin_t,
              Vector& surface_scalar_reflectivity,
              ConstVectorView scat_za_grid,
              const Index& nstreams,
              const Index& non_iso_inc,
              const String& pfct_method,
              const Verbosity& verbosity )
{
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
  Vector scat_angle_grid;
  Index pfct_za_grid_size;
  if( pfct_method=="interpolate" )
  {
    pfct_za_grid_size=181;
    nlinspace(scat_angle_grid, 0, 180, pfct_za_grid_size);
  }
  else
    // Scattering angle grid, assumed here that it is the same for
    // all scattering elements
    scat_angle_grid = scat_data[0][0].za_grid;
  Matrix phase_function(nlyr,scat_angle_grid.nelem(), 0.);
  
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
  for (Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
    {
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

//#pragma omp critical(fortran_disort)
//      {
      dtauc_ssalbCalc(ws, dtauc, ssalb, scat_data, f_index,
                      propmat_clearsky_agenda,
                      pnd_field,
                      t_field(Range(0,nlyr+1),joker,joker),
                      z_field(Range(0,nlyr+1),joker,joker),
                      vmr_field(joker,Range(0,nlyr+1),joker,joker),
                      p_grid[Range(0,nlyr+1)],
                      cloudbox_limits, f_grid[Range(f_index,1)], verbosity);

      if( pfct_method=="interpolate" )
      {
        phase_functionCalc2(phase_function,
                            scat_data, f_index,
                            pnd_field, t_field, cloudbox_limits,
                            pfct_za_grid_size, verbosity);
        for( Index l=0; l<nlyr; l++ )
          if( phase_function(l,0)==0. )
            assert( ssalb[l]==0. );

        //cout << "entering pmomCalc for f_index=" << f_index << " (f="
        //     << f_grid[f_index]*1e-9 << "GHz).\n";
        pmomCalc2(pmom, phase_function, scat_angle_grid, n_legendre, verbosity);
      }
      else
      {
        phase_functionCalc(phase_function, scat_data, f_index, pnd_field,
                           cloudbox_limits, pfct_method );
        for( Index l=0; l<nlyr; l++ )
          if( phase_function(l,0)==0. )
            assert( ssalb[l]==0. );

        //cout << "entering pmomCalc for f_index=" << f_index << " (f="
        //     << f_grid[f_index]*1e-9 << "GHz).\n";
        pmomCalc(pmom, phase_function, scat_angle_grid, n_legendre, verbosity);
      }

      // Wavenumber in [1/cm]
      Numeric wvnmlo = f_grid[f_index]/(100*SPEED_OF_LIGHT);
      Numeric wvnmhi = wvnmlo;
      
      // calculate radiant quantities at boundary of computational layers. 
      Index usrtau = FALSE_; 
      
      //DEBUG_VAR(dtauc)
      
      //cout << "entering Disort calc at freq[" << f_index << "]="
      //     << f_grid[f_index]*1e-9 << " GHz\n"
      //     << "  with surftemp=" << surface_skin_t
      //     << " K and albedo=" << surface_scalar_reflectivity[f_index]
      //     << "\n";
      
// JM: once (2-3-454), I extended the critical region due to
// modified-variable-issues inside Disort. However, later on I couldn't
// reproduce the problems anymore. Since the extension created problems in catching
// errors thrown inside methods called within the critical region, the region is
// reduced to its original extend (2-3-486), covering only the Disort call
// itself. If any kind of fishy behaviour is observed, we have to reconsider
// extending (and proper error handling) again.
#pragma omp critical(fortran_disort)
      {
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
                  &surface_scalar_reflectivity[f_index], hl.get_c_array(),
                  &surface_skin_t, &ttemp, &temis,
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
      }

      for(Index j = 0; j<numu; j++)
          for(Index k = 0; k<(cloudbox_limits[1]-cloudbox_limits[0]+1); k++)
            doit_i_field(f_index, k, 0, 0, j, 0, 0) =
              uu(0,nlyr-k-cloudbox_limits[0],j) / (100*SPEED_OF_LIGHT);
    }
  delete [] prnt;
}


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
#pragma omp critical(fortran_disort_qgausn)
  {
      qgausn_(&nn,
              gmu.get_c_array(),
              gwt.get_c_array()
              );
  }

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

void run_disort( Workspace&,
              Tensor7&,
              ConstVectorView,
              ConstVectorView,
              ConstTensor3View,
              ConstTensor3View,
              ConstTensor4View,
              ConstTensor4View,
              const ArrayOfArrayOfSingleScatteringData&,
              const Agenda&,
              const Agenda&,
              const ArrayOfIndex&,
              Numeric&,
              Vector&,
              ConstVectorView,
              const Index&,
              const Index&,
              const String&,
              const Verbosity& )
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
}

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


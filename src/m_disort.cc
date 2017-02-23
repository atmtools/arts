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
#include "math_funcs.h"
#include "messages.h"
#include "m_general.h"
#include "wsv_aux.h"
#include "xml_io.h"


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
                const Numeric& surface_skin_t,
                const Vector& surface_scalar_reflectivity,
                const Index& nstreams,
                const Index& non_iso_inc,
                const String& pfct_method,
                const Verbosity& verbosity )
{
  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  // FIXME: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  check_disort_input( cloudbox_on, disort_is_initialized,
                      atmfields_checked, atmgeom_checked, cloudbox_checked,
                      scat_data, scat_za_grid, nstreams, pfct_method,
                      pnd_field.ncols(), doit_i_field.npages() );

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  get_disortsurf_props( albedo, btemp,
                        f_grid, surface_skin_t, surface_scalar_reflectivity );

  run_disort( ws, doit_i_field,
              f_index, f_grid,p_grid,z_field, t_field,vmr_field, pnd_field,
              scat_data, scat_data_mono,
              propmat_clearsky_agenda, opt_prop_part_agenda, spt_calc_agenda,
              iy_main_agenda,
              cloudbox_limits,
              btemp, albedo,
              scat_za_grid, nstreams,
              non_iso_inc, pfct_method,
              verbosity );
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcWithARTSSurface(Workspace& ws,
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
                const Agenda& surface_rtprop_agenda,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Vector& scat_za_grid,
                const Index& nstreams,
                const Index& non_iso_inc,
                const String& pfct_method,
                const Verbosity& verbosity )
{
  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  // FIXME: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  check_disort_input( cloudbox_on, disort_is_initialized,
                      atmfields_checked, atmgeom_checked, cloudbox_checked,
                      scat_data, scat_za_grid, nstreams, pfct_method,
                      pnd_field.ncols(), doit_i_field.npages() );

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  // for now, surface at lowest atm level. later use z_surface or the like
  // for that.
  // at the moment this is only required for groundtype "A", but 
  const Numeric surf_altitude = z_field(0,0,0);
  //const Numeric surf_altitude = z_surface(0,0);

  surf_albedoCalc( ws, albedo, btemp,
                   surface_rtprop_agenda,
                   f_grid, scat_za_grid, surf_altitude,
                   verbosity );

  run_disort( ws, doit_i_field,
              f_index, f_grid,p_grid,z_field, t_field,vmr_field, pnd_field,
              scat_data, scat_data_mono,
              propmat_clearsky_agenda, opt_prop_part_agenda, spt_calc_agenda,
              iy_main_agenda,
              cloudbox_limits,
              btemp, albedo,
              scat_za_grid, nstreams,
              non_iso_inc, pfct_method,
              verbosity );
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
                const Numeric&,
                const Vector&,
                const Index&,
                const Index&,
                const String&,
                const Verbosity&)
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
}

void DisortCalcWithARTSSurface(Workspace&,
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
                const Agenda&,
                const Tensor4&,
                const Tensor3&,
                const Tensor3&,
                const Tensor4&,
                const Vector&,
                const ArrayOfArrayOfSingleScatteringData&,
                const Vector&,
                const Vector&,
                const Index&,
                const Index&,
                const String&,
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
              const Verbosity& verbosity _U_ )
{
  if (!cloudbox_on)
  {
    //CREATE_OUT0;
    //disort_is_initialized = 0;
    //out0 << "  Cloudbox is off, scattering calculation will be skipped.\n";
    //return;
    throw runtime_error( "Cloudbox is off, no scattering calculations to be"
                         "performed." );
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
      os << "DISORT calculations currently only possible with "
         << "lower cloudbox limit\n"
         << "at 0th atmospheric level "
         << "(assumes surface there, ignoring z_surface).\n";
      throw runtime_error(os.str());
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
    
  //------------- end of checks ---------------------------------------
  
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index Nza = scat_za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  doit_i_field.resize( Nf, Np_cloud, 1, 1, Nza, 1, stokes_dim );
  doit_i_field = NAN;
  
  disort_is_initialized = 1;
}


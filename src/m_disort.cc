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
                // WS Input
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& scat_data_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits,
                const Agenda& propmat_clearsky_agenda, 
                const Index& atmosphere_dim,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Vector& scat_za_grid,
                const Index& stokes_dim,
                const Numeric& surface_skin_t,
                const Vector& surface_scalar_reflectivity,
                const Index& nstreams,
                const String& pfct_method,
                const Index& new_optprop,
                const Index& Npfct,
                const Verbosity& verbosity )
{
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on)
  {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  // FIXME: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  check_disort_input( cloudbox_on,
                      atmfields_checked, atmgeom_checked,
                      cloudbox_checked, scat_data_checked,
                      atmosphere_dim, stokes_dim, cloudbox_limits,
                      scat_data, scat_za_grid, nstreams, pfct_method,
                      pnd_field.ncols() );

  init_ifield( doit_i_field, f_grid, cloudbox_limits, scat_za_grid.nelem(), stokes_dim );

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;


  get_disortsurf_props( albedo, btemp,
                        f_grid, surface_skin_t, surface_scalar_reflectivity );

  if( new_optprop )
    run_disort2( ws, doit_i_field,
              f_grid,
              p_grid, z_field, t_field, vmr_field, pnd_field,
              scat_data,
              propmat_clearsky_agenda,
              cloudbox_limits,
              btemp, albedo,
              scat_za_grid, nstreams,
              pfct_method, Npfct,
              verbosity );
  else
    run_disort( ws, doit_i_field,
              f_grid,
              p_grid, z_field, t_field, vmr_field, pnd_field,
              scat_data,
              propmat_clearsky_agenda,
              cloudbox_limits,
              btemp, albedo,
              scat_za_grid, nstreams,
              pfct_method,
              verbosity );
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcWithARTSSurface(Workspace& ws,
                // WS Output:
                Tensor7& doit_i_field,
                // WS Input
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& scat_data_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits, 
                const Agenda& propmat_clearsky_agenda, 
                const Agenda& surface_rtprop_agenda,
                const Index& atmosphere_dim,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Vector& scat_za_grid,
                const Index& stokes_dim,
                const Index& nstreams,
                const String& pfct_method,
                const Index& new_optprop,
                const Index& Npfct,
                const Verbosity& verbosity )
{
  if (!cloudbox_on)
  {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  // FIXME: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  check_disort_input( cloudbox_on,
                      atmfields_checked, atmgeom_checked,
                      cloudbox_checked, scat_data_checked,
                      atmosphere_dim, stokes_dim, cloudbox_limits,
                      scat_data, scat_za_grid, nstreams, pfct_method,
                      pnd_field.ncols() );

  init_ifield( doit_i_field, f_grid, cloudbox_limits, scat_za_grid.nelem(), stokes_dim );

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

  if( new_optprop )
      run_disort2( ws, doit_i_field,
                   f_grid,
                   p_grid, z_field, t_field, vmr_field, pnd_field,
                   scat_data,
                   propmat_clearsky_agenda,
                   cloudbox_limits,
                   btemp, albedo,
                   scat_za_grid, nstreams,
                   pfct_method, Npfct,
                   verbosity );
  else
      run_disort( ws, doit_i_field,
                  f_grid,
                  p_grid, z_field, t_field, vmr_field, pnd_field,
                  scat_data,
                  propmat_clearsky_agenda,
                  cloudbox_limits,
                  btemp, albedo,
                  scat_za_grid, nstreams,
                  pfct_method,
                  verbosity );
}

#else /* ENABLE_DISORT */

void DisortCalc(Workspace&,
                // WS Output:
                Tensor7&,
                // WS Input
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const ArrayOfIndex&,
                const Agenda&,
                const Index&,
                const Tensor4&,
                const Tensor3&,
                const Tensor3&,
                const Tensor4&,
                const Vector&,
                const ArrayOfArrayOfSingleScatteringData&,
                const Vector&,
                const Vector&,
                const Index&,
                const Numeric&,
                const Vector&,
                const Index&,
                const String&,
                const Index&,
                const Index&,
                const Verbosity&)
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
}

void DisortCalcWithARTSSurface(Workspace&,
                // WS Output:
                Tensor7&,
                // WS Input
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const ArrayOfIndex&,
                const Agenda&,
                const Agenda&,
                const Index&,
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
                const Index&,
                const Index&,
                const Verbosity&)
{
  throw runtime_error ("This version of ARTS was compiled without DISORT support.");
}

#endif /* ENABLE_DISORT */


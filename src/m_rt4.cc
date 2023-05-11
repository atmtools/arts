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

#include <complex.h>
#include <stdexcept>
#include "auto_md.h"
#include "disort.h"
#include "m_xml.h"
#include "rt4.h"
#include "species_tags.h"

#ifdef ENABLE_RT4

const Numeric pfct_threshold = 0.05;

/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Calc(Workspace& ws,
             // WS Output:
             Tensor7& cloudbox_field,
             Vector& za_grid,
             Vector& aa_grid,
             // WS Input
             const Index& atmfields_checked,
             const Index& atmgeom_checked,
             const Index& scat_data_checked,
             const Index& cloudbox_checked,
             const Index& cloudbox_on,
             const ArrayOfIndex& cloudbox_limits,
             const Agenda& propmat_clearsky_agenda,
             const Agenda& surface_rtprop_agenda,
             const Tensor4& pnd_field,
             const AtmField& atm_field,
             const ArrayOfArrayOfSingleScatteringData& scat_data,
             const ArrayOfArrayOfSpeciesTag& abs_species,
             const Vector& f_grid,
             const Index& stokes_dim,
             const SurfaceField& surface_field,
             const Index& nstreams,
             const String& pfct_method,
             const String& quadtype,
             const Index& add_straight_angles,
             const Index& pfct_aa_grid_size,
             //const Numeric& pfct_threshold,
             const Index& auto_inc_nstreams,
             const Index& robust,
             const Index& za_interp_order,
             const Index& cos_za_interp,
             const Numeric& max_delta_tau) {
  if (!cloudbox_on) {
    return;
  }

  //! Fixme: Must have lat and lon
  const Numeric surf_alt = surface_field.single_value(Surf::Key::h, 0, 0);

  const String quad_type = quadtype.toupper();
  Index nhza, nhstreams, nummu;
  
  check_rt4_input(nhstreams,
                  nhza,
                  nummu,
                  cloudbox_on,
                  atmfields_checked,
                  atmgeom_checked,
                  cloudbox_checked,
                  scat_data_checked,
                  cloudbox_limits,
                  scat_data,
                  stokes_dim,
                  nstreams,
                  quad_type,
                  add_straight_angles,
                  pnd_field.ncols());

  init_ifield(cloudbox_field, f_grid, cloudbox_limits, 2 * nummu, 1, stokes_dim);

  // in RT4 mu_values is generally only output. however, we need the values for
  // preparing the single scattering data at these angles. therefore, we
  // calculate them here using RT4's proprietary quadrature methods. They
  // simultaneously provide the quadrature weights, too. We keep them so far,
  // might use them for ensuring proper normalization in the preparation of the
  // single scattering data.
  Vector mu_values(nummu, 0.);
  Vector quad_weights(nummu, 0.);

  get_quad_angles(mu_values,
                  quad_weights,
                  za_grid,
                  aa_grid,
                  quad_type,
                  nhstreams,
                  nhza,
                  nummu);

  // Preparing surface setup.
  //

  // Initializing surface related interface-related RT4 interface parameters.
  const Index nf = f_grid.nelem();

  // dummy values for parameters not relevant for this ground_type
  Numeric surface_skin_t = 0.;
  Vector ground_albedo(nf, 0.);
  Tensor3 ground_reflec(nf, stokes_dim, stokes_dim, 0.);
  Complex gidef(1, 0.);
  ComplexVector ground_index(nf, gidef);

  // parameters that will be updated below
  Tensor5 surf_refl_mat(nf, nummu, stokes_dim, nummu, stokes_dim, 0.);
  Tensor3 surf_emis_vec(nf, nummu, stokes_dim, 0.);

  surf_optpropCalc(ws,
                   surf_refl_mat,
                   surf_emis_vec,
                   surface_rtprop_agenda,
                   f_grid,
                   za_grid,
                   mu_values,
                   quad_weights,
                   stokes_dim,
                   surf_alt);

  run_rt4(ws,
          cloudbox_field,
          za_grid,
          f_grid,
          atm_field,
          pnd_field(joker, joker, 0, 0),
          scat_data,
          abs_species,
          propmat_clearsky_agenda,
          cloudbox_limits,
          stokes_dim,
          nummu,
          nhza,
          "A",
          surface_skin_t,
          ground_albedo,
          ground_reflec,
          ground_index,
          surf_refl_mat,
          surf_emis_vec,
          surface_rtprop_agenda,
          surf_alt,
          quad_type,
          mu_values,
          quad_weights,
          auto_inc_nstreams,
          robust,
          za_interp_order,
          cos_za_interp,
          pfct_method,
          pfct_aa_grid_size,
          pfct_threshold,
          max_delta_tau);

  za_grid_adjust(za_grid, mu_values, nummu);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RT4CalcWithRT4Surface(Workspace& ws,
                           // WS Output:
                           Tensor7& cloudbox_field,
                           Vector& za_grid,
                           Vector& aa_grid,
                           // WS Input
                           const Index& atmfields_checked,
                           const Index& atmgeom_checked,
                           const Index& scat_data_checked,
                           const Index& cloudbox_checked,
                           const Index& cloudbox_on,
                           const ArrayOfIndex& cloudbox_limits,
                           const Agenda& propmat_clearsky_agenda,
                           const Tensor4& pnd_field,
                           const AtmField& atm_field,
                           const ArrayOfArrayOfSingleScatteringData& scat_data,
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const Vector& f_grid,
                           const Index& stokes_dim,
                           const SurfaceField& surface_field,
                           const Numeric& surface_skin_t,
                           const Vector& surface_scalar_reflectivity,
                           const Tensor3& surface_reflectivity,
                           const GriddedField3& surface_complex_refr_index,
                           const Index& nstreams,
                           const String& pfct_method,
                           const String& groundtype,
                           const String& quadtype,
                           const Index& add_straight_angles,
                           const Index& pfct_aa_grid_size,
                           //const Numeric& pfct_threshold,
                           const Index& auto_inc_nstreams,
                           const Index& robust,
                           const Index& za_interp_order,
                           const Index& cos_za_interp,
                           const Numeric& max_delta_tau) {
  if (!cloudbox_on) {
    return;
  }

  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  const String quad_type = quadtype.toupper();
  Index nhza, nhstreams, nummu;
  check_rt4_input(nhstreams,
                  nhza,
                  nummu,
                  cloudbox_on,
                  atmfields_checked,
                  atmgeom_checked,
                  cloudbox_checked,
                  scat_data_checked,
                  cloudbox_limits,
                  scat_data,
                  stokes_dim,
                  nstreams,
                  quad_type,
                  add_straight_angles,
                  pnd_field.ncols());

  init_ifield(cloudbox_field, f_grid, cloudbox_limits, 2 * nummu, 1, stokes_dim);

  // in RT4 mu_values is generally only output. however, we need the values for
  // preparing the single scattering data at these angles. therefore, we
  // calculate them here using RT4's proprietary quadrature methods. They
  // simultaneously provide the quadrature weights, too. We keep them so far,
  // might use them for ensuring proper normalization in the preparation of the
  // single scattering data.
  Vector mu_values(nummu, 0.);
  Vector quad_weights(nummu, 0.);

  get_quad_angles(mu_values,
                  quad_weights,
                  za_grid,
                  aa_grid,
                  quad_type,
                  nhstreams,
                  nhza,
                  nummu);

  // Preparing surface setup.
  //

  // Initializing surface related interface-related RT4 interface parameters.
  const Index nf = f_grid.nelem();
  const String ground_type = groundtype.toupper();

  // dummy values for parameters not relevant for this ground_type
  Tensor5 surf_refl_mat(nf, nummu, stokes_dim, nummu, stokes_dim, 0.);
  Tensor3 surf_emis_vec(nf, nummu, stokes_dim, 0.);

  // parameters that will be updated below
  Vector ground_albedo(nf, 0.);
  Tensor3 ground_reflec(nf, stokes_dim, stokes_dim, 0.);
  Complex gidef(1, 0.);
  ComplexVector ground_index(nf, gidef);

  get_rt4surf_props(ground_albedo,
                    ground_reflec,
                    ground_index,
                    f_grid,
                    ground_type,
                    surface_skin_t,
                    surface_scalar_reflectivity,
                    surface_reflectivity,
                    surface_complex_refr_index,
                    stokes_dim);

  Agenda dummy_agenda{ws};

  run_rt4(ws,
          cloudbox_field,
          za_grid,
          f_grid,
          atm_field,
          pnd_field(joker, joker, 0, 0),
          scat_data,
          abs_species,
          propmat_clearsky_agenda,
          cloudbox_limits,
          stokes_dim,
          nummu,
          nhza,
          ground_type,
          surface_skin_t,
          ground_albedo,
          ground_reflec,
          ground_index,
          surf_refl_mat,
          surf_emis_vec,
          dummy_agenda,
          surface_field.single_value(Surf::Key::h, 0, 0),  // FIXME: Must know lat and lon
          quad_type,
          mu_values,
          quad_weights,
          auto_inc_nstreams,
          robust,
          za_interp_order,
          cos_za_interp,
          pfct_method,
          pfct_aa_grid_size,
          pfct_threshold,
          max_delta_tau);

  za_grid_adjust(za_grid, mu_values, nummu);
}

#else /* ENABLE_RT4 */
/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Calc(Workspace& ws,
             // WS Output:
             Tensor7& cloudbox_field,
             Vector& za_grid,
             Vector& aa_grid,
             // WS Input
             const Index& atmfields_checked,
             const Index& atmgeom_checked,
             const Index& scat_data_checked,
             const Index& cloudbox_checked,
             const Index& cloudbox_on,
             const ArrayOfIndex& cloudbox_limits,
             const Agenda& propmat_clearsky_agenda,
             const Agenda& surface_rtprop_agenda,
             const Tensor4& pnd_field,
             const AtmField& atm_field,
             const ArrayOfArrayOfSingleScatteringData& scat_data,
             const ArrayOfArrayOfSpeciesTag& abs_species,
             const Vector& f_grid,
             const Index& stokes_dim,
             const SurfaceField& surface_field,
             const Index& nstreams,
             const String& pfct_method,
             const String& quadtype,
             const Index& add_straight_angles,
             const Index& pfct_aa_grid_size,
             //const Numeric& pfct_threshold,
             const Index& auto_inc_nstreams,
             const Index& robust,
             const Index& za_interp_order,
             const Index& cos_za_interp,
             const Numeric& max_delta_tau) {
  ARTS_USER_ERROR("This version of ARTS was compiled without RT4 support.");
}
/* Workspace method: Doxygen documentation will be auto-generated */
void RT4CalcWithRT4Surface(Workspace& ws,
                           // WS Output:
                           Tensor7& cloudbox_field,
                           Vector& za_grid,
                           Vector& aa_grid,
                           // WS Input
                           const Index& atmfields_checked,
                           const Index& atmgeom_checked,
                           const Index& scat_data_checked,
                           const Index& cloudbox_checked,
                           const Index& cloudbox_on,
                           const ArrayOfIndex& cloudbox_limits,
                           const Agenda& propmat_clearsky_agenda,
                           const Tensor4& pnd_field,
                           const AtmField& atm_field,
                           const ArrayOfArrayOfSingleScatteringData& scat_data,
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const Vector& f_grid,
                           const Index& stokes_dim,
                           const SurfaceField& surface_field,
                           const Numeric& surface_skin_t,
                           const Vector& surface_scalar_reflectivity,
                           const Tensor3& surface_reflectivity,
                           const GriddedField3& surface_complex_refr_index,
                           const Index& nstreams,
                           const String& pfct_method,
                           const String& groundtype,
                           const String& quadtype,
                           const Index& add_straight_angles,
                           const Index& pfct_aa_grid_size,
                           //const Numeric& pfct_threshold,
                           const Index& auto_inc_nstreams,
                           const Index& robust,
                           const Index& za_interp_order,
                           const Index& cos_za_interp,
                           const Numeric& max_delta_tau) {
  ARTS_USER_ERROR("This version of ARTS was compiled without RT4 support.");
}

#endif /* ENABLE_RT4 */

/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef ENABLE_RT4
void RT4Test(Tensor4& out_rad,
             const String& datapath) {
  rt4_test(out_rad, datapath);
}
#else
void RT4Test(Tensor4&, const String&) {
  ARTS_USER_ERROR ("This version of ARTS was compiled without RT4 support.");
}
#endif

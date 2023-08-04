/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   sun.cc
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Functions needed for radiative transfer with direct sources.
*/

#include "sun.h"
#include "atm.h"
#include <workspace.h>
#include "arts_conversions.h"
#include "check_input.h"
#include "debug.h"
#include "matpack_data.h"
#include "physics_funcs.h"
#include "geodetic.h"
#include "arts.h"
#include "surf.h"

using Constant::sigma;
using Constant::pi;

/*===========================================================================
  === The functions
  ===========================================================================*/

std::ostream& operator<<(std::ostream& os, const Sun& sun) {
  os << "Sun: " << sun.description;
  os << " Radius: " << sun.radius << "m ";
  os << " Distance: " << sun.distance << "m \n";
  os << " Latitude: " << sun.latitude << "° \n";
  os << " Longitude: " << sun.longitude << "° \n";
  os << " Spectrum [W/m2/Hz]: \n" << sun.spectrum ;
  return os;
}

void get_scattered_sunsource(Workspace& ws,
                              StokvecVector& scattered_sunlight,
                              const Vector& f_grid,
                              const AtmPoint& atm_point,
                              const Matrix& transmitted_sunlight,
                              const Vector& gas_scattering_los_in,
                              const Vector& gas_scattering_los_out,
                              const Agenda& gas_scattering_agenda) {
  PropmatVector K_sca;
  MuelmatVector gas_scattering_mat;
  Vector sca_fct_dummy;

  // calculate gas scattering properties
  gas_scattering_agendaExecute(ws,
                               K_sca,
                               gas_scattering_mat,
                               sca_fct_dummy,
                               f_grid,
                               atm_point,
                               gas_scattering_los_in,
                               gas_scattering_los_out,
                               0,
                               gas_scattering_agenda);

  //some basic quantities
  Index ns = transmitted_sunlight.ncols();
  Index nf = f_grid.nelem();

  //allocate and resize
  StokvecVector scattered_sunlight_temp(1, ns);

  Matrix mat_temp(1, ns,0.);
  // Calculate the scattered radiation
  for (Index i_f = 0; i_f < nf; i_f++) {
    Stokvec scattered_sunlight_temp{transmitted_sunlight[i_f]};
    scattered_sunlight[i_f] = K_sca[i_f].A() / (4*pi) * gas_scattering_mat[i_f] * scattered_sunlight_temp;
  }

  //TODO: Include jacobian mechanism
}

void get_sun_background(Matrix& iy,
                         Index& suns_visible,
                         const ArrayOfSun& suns,
                         const Ppath& ppath,
                         const Vector& f_grid,
                         const Vector2 refellipsoid) {
  const Index np = ppath.np;

  //set visibilty flag to default
  suns_visible = 0;

  //allocate iy and set it to zero
  iy.resize(f_grid.nelem(), 4);
  iy=0.;

  Vector rtp_pos, rtp_los;
  rtp_pos.resize(3);
  rtp_pos = ppath.pos(np - 1, Range(0, 3));
  rtp_los.resize(ppath.los.ncols());
  rtp_los = ppath.los(np - 1, joker);

  // if (ppath_what_background(ppath) == 9 || ppath_what_background(ppath) == 1){
  //   for (Index i_sun = 0; i_sun < suns.nelem(); i_sun++) {
  //     get_sun_radiation(iy, suns_visible, suns[i_sun], rtp_pos, rtp_los, refellipsoid);
  //   }
  // }
  ARTS_USER_ERROR("ERROR")
}

void get_sun_radiation(Matrix& iy,
                        Index& suns_visible,
                         const Sun& sun,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector2 refellipsoid) {

  //Calculate earth centric radial component of sun_pos and rtp_pos.
  const Numeric R_sun = sun.distance;
  const Numeric R_rte = 0;ARTS_USER_ERROR("ERROR")//refell2r(refellipsoid, rtp_pos[1]) + rtp_pos[0];

  //Transform to cartesian coordinate system
  Numeric r_sun_x, r_sun_y, r_sun_z;
  Numeric r_rte_x, r_rte_y, r_rte_z;
  Numeric r_los_x, r_los_y, r_los_z;

  // r_sun
  //sph2cart(r_sun_x, r_sun_y, r_sun_z, R_sun, sun.latitude, sun.longitude);
ARTS_USER_ERROR("ERROR")
  // r_rte, r_los
  // poslos2cart(r_rte_x,
  //             r_rte_y,
  //             r_rte_z,
  //             r_los_x,
  //             r_los_y,
  //             r_los_z,
  //             R_rte,
  //             rtp_pos[1],
  //             rtp_pos[2],
  //             rtp_los[0],
  //             rtp_los[1]);
  ARTS_USER_ERROR("ERROR")

  //Calculate vector of line of sight and unit vector pointing from
  //ppath point to the sun.
  const Numeric r_ps_x = r_sun_x - r_rte_x;
  const Numeric r_ps_y = r_sun_y - r_rte_y;
  const Numeric r_ps_z = r_sun_z - r_rte_z;

  //abs value of r_ps
  const Numeric r_ps =
      sqrt(r_ps_x * r_ps_x + r_ps_y * r_ps_y + r_ps_z * r_ps_z);

  //abs value of r_los
  const Numeric r_glos =
      sqrt(r_los_x * r_los_x + r_los_y * r_los_y + r_los_z * r_los_z);

  //Calculate angle beta between line of sight and the line between ppath point and the sun
  //using scalar product
  const Numeric cos_beta =
      (r_ps_x * r_los_x + r_ps_y * r_los_y + r_ps_z * r_los_z) / (r_ps * r_glos);
  const Numeric beta = acos(cos_beta);

  // angular radius of sun
  const Numeric alpha = atan(sun.radius / r_ps);

  //Check if we see the sun. We see the sun if the angle beta is smaller than
  // the angular radius alpha of the sun.
  if (beta <= alpha) {
    //Here we assume that the sun radiates isotropically.
    Matrix sun_radiance = sun.spectrum;
    sun_radiance /= pi;

    iy += sun_radiance;

    // visibility flag
    suns_visible = 1;
  }
}

void get_direct_radiation(Workspace& ws,
                          ArrayOfMatrix& direct_radiation,
                          ArrayOfArrayOfTensor3& ddirect_radiation_dx,
                          const Vector& f_grid,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const AtmField& atm_field,
                          const Index& cloudbox_on,
                          const ArrayOfIndex& cloudbox_limits,
                          const Index& gas_scattering_do,
                          const Index& irradiance_flag,
                          const ArrayOfPpath& sun_ppaths,
                          const ArrayOfSun& suns,
                          const ArrayOfIndex& suns_visible,
                          const Vector2 refellipsoid,
                          const Tensor4& pnd_field,
                          const ArrayOfTensor4& dpnd_field_dx,
                          const ArrayOfString& scat_species,
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const Index& jacobian_do,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Agenda& propmat_clearsky_agenda,
                          const Agenda& water_p_eq_agenda,
                          const Agenda& gas_scattering_agenda,
                          const Numeric& rte_alonglos_v) {
  Vector sun_pos(3);

  //dummy variables needed for the output and input of iyTransmission and
  // gas_scattering_agenda
  ArrayOfMatrix iy_aux_dummy;
  ArrayOfString iy_aux_vars_dummy;
  ArrayOfAtmPoint ppvar_atm_dummy;
  Matrix ppvar_pnd_dummy;
  ArrayOfVector ppvar_f_dummy;
  Tensor3 ppvar_iy_dummy;
  Tensor4 ppvar_trans_cumulat_dummy;
  Tensor4 ppvar_trans_partial_dummy;

  direct_radiation.resize(suns.nelem(),Matrix(f_grid.nelem(), 4, 0.));

  Matrix radiation_toa;
  Matrix radiation_trans;
  ArrayOfTensor3 dradiation_trans;
  Vector rtp_pos, rtp_los;
  Index np;

ARTS_USER_ERROR("ERROR")
//  const auto& lat_grid = atm_field.grid[1];
  //const auto& lon_grid = atm_field.grid[2];
Vector lat_grid, lon_grid;

  for (Index i_sun = 0; i_sun < suns.nelem(); i_sun++) {
    np = sun_ppaths[i_sun].np;

    if (suns_visible[i_sun]) {
      sun_pos = {suns[i_sun].distance,
                  suns[i_sun].latitude,
                  suns[i_sun].longitude};

      // we need the distance to the sun relative to the surface
      // sun_pos[0] =
      //     sun_pos[0] -
      //     pos2refell_r(
      //         3, refellipsoid, lat_grid, lon_grid, sun_pos);
ARTS_USER_ERROR("ERROR")

      if (irradiance_flag) {
        //Get the spectral irradiance

        //get the TOA distance to earth center.
        Numeric R_TOA;// =
            // refell2r(refellipsoid, sun_ppaths[i_sun].pos(np - 1, 1)) +
            // sun_ppaths[i_sun].pos(np - 1, 0);
ARTS_USER_ERROR("ERROR")
        //get the distance between sun and sun_ppath at TOA
        Numeric R_Sun2Toa;
        // distance3D(R_Sun2Toa,
        //            R_TOA,
        //            sun_ppaths[i_sun].pos(np - 1, 1),
        //            sun_ppaths[i_sun].pos(np - 1, 2),
        //            sun_pos[0],
        //            sun_pos[1],
        //            sun_pos[2]);
ARTS_USER_ERROR("ERROR")

        // Scale the incoming sun_irradiance spectrum
        radiation_toa = suns[i_sun].spectrum;
        radiation_toa *= suns[i_sun].radius * suns[i_sun].radius;
        radiation_toa /= (suns[i_sun].radius * suns[i_sun].radius +
                      R_Sun2Toa * R_Sun2Toa);
      } else {
        //Get the spectral radiance instead

        //Set sun position
        rtp_pos.resize(3);
        rtp_pos = sun_ppaths[i_sun].pos(np - 1, Range(0, 3));
        rtp_los.resize(sun_ppaths[i_sun].los.ncols());
        rtp_los = sun_ppaths[i_sun].los(np - 1, joker);

        radiation_toa.resize(f_grid.nelem(), 4);
        radiation_toa = 0;

        Index visible;
        get_sun_radiation(radiation_toa, visible, suns[i_sun], rtp_pos, rtp_los, refellipsoid);
      }


      iyTransmissionStandard(ws,
                             radiation_trans,
                             iy_aux_dummy,
                             dradiation_trans,
                             ppvar_atm_dummy,
                             ppvar_pnd_dummy,
                             ppvar_f_dummy,
                             ppvar_iy_dummy,
                             ppvar_trans_cumulat_dummy,
                             ppvar_trans_partial_dummy,
                             f_grid,
                             abs_species,
                             atm_field,
                             cloudbox_on,
                             cloudbox_limits,
                             gas_scattering_do,
                             pnd_field,
                             dpnd_field_dx,
                             scat_species,
                             scat_data,
                             iy_aux_vars_dummy,
                             jacobian_do,
                             jacobian_quantities,
                             sun_ppaths[i_sun],
                             radiation_toa,
                             propmat_clearsky_agenda,
                             water_p_eq_agenda,
                             gas_scattering_agenda,
                             1,
                             Tensor3(),
                             rte_alonglos_v);

      if (jacobian_do && dradiation_trans.nelem()){
        ddirect_radiation_dx[i_sun] = dradiation_trans;
      }
      direct_radiation[i_sun] = radiation_trans;
    }
  }
}

void get_sun_ppaths(Workspace& ws,
                     ArrayOfPpath& sun_ppaths,
                     ArrayOfIndex& suns_visible,
                     ArrayOfVector& sun_rte_los,
                     const Vector& rte_pos,
                     const ArrayOfSun& suns,
                     const Vector& f_grid,
                     const AtmField& atm_field,
                     const SurfaceField& surface_field,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Agenda& ppath_step_agenda) {
ARTS_USER_ERROR("ERROR")
//const auto& lat_grid = atm_field.grid[1];
//const auto& lon_grid = atm_field.grid[2];
Vector lat_grid, lon_grid;

  Vector sun_pos(3);
  Vector sun_rte_los_isun(2);
  Numeric ppath_lraytrace2 = ppath_lraytrace;
  Ppath sun_ppath;

  for (Index i_sun = 0; i_sun < suns.nelem(); i_sun++) {
    sun_pos = {suns[i_sun].distance,
                suns[i_sun].latitude,
                suns[i_sun].longitude};

    // we need the distance to the sun relative to the surface
    // sun_pos[0] =
    //     sun_pos[0] -
    //     pos2refell_r(
    //         3, refellipsoid, lat_grid, lon_grid, sun_pos);
ARTS_USER_ERROR("ERROR")

    // get the line of sight direction from sun i_sun to ppath point
    ARTS_ASSERT(false)
    /*rte_losGeometricFromRtePosToRtePos2(sun_rte_los_isun,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        rte_pos,
                                        sun_pos,
                                        );*/

    // calculate ppath (sun path) from sun to ppath point
    ARTS_ASSERT(false)
    /*FIXME
    ppathFromRtePos2(ws,
                     sun_ppath,
                     sun_rte_los_isun,
                     ppath_lraytrace2,
                     ppath_step_agenda,
                     atm_field,
                     f_grid,
                     refellipsoid,
                     z_surface,
                     rte_pos,
                     sun_pos,
                     ppath_lmax,
                     2e-5,
                     5.,
                     0.5,
                     );
                     */

    sun_ppaths[i_sun] = sun_ppath;
    sun_rte_los[i_sun] = sun_rte_los_isun;

    // Check that sun path hast more than 1 ppath points and that space is background
    // suns_visible[i_sun] =
    //     sun_ppath.np > 1 && ppath_what_background(sun_ppath) == 9;
    ARTS_USER_ERROR("ERROR")
  }
}

Matrix regrid_sun_spectrum(const GriddedField2 &sun_spectrum_raw,
                          const Vector &f_grid,
                          const Numeric &temperature){
  const Index nf = f_grid.nelem();
  const Vector data_f_grid = sun_spectrum_raw.get_numeric_grid(0);
  const Numeric data_fmin = data_f_grid[0];
  const Numeric data_fmax = data_f_grid[data_f_grid.nelem() - 1];

    // Result array
  Matrix int_data(f_grid.nelem(), 4, 0.);

  const Numeric* f_grid_begin = f_grid.unsafe_data_handle();
  const Numeric* f_grid_end = f_grid_begin + f_grid.nelem();
  const Index i_fstart = std::distance(
      f_grid_begin, std::lower_bound(f_grid_begin, f_grid_end, data_fmin));
  const Index i_fstop =
      std::distance(
          f_grid_begin,
          std::upper_bound(f_grid_begin + i_fstart, f_grid_end, data_fmax)) -
      1;

  // Ignore band if all frequencies are below or above data_f_grid:
  if (i_fstart == nf || i_fstop == -1) {
  }

  const Index f_extent = i_fstop - i_fstart + 1;

  // If f_extent is less than one, then the entire data_f_grid is between two
  // grid points of f_grid. (So that we do not have any f_grid points inside
  // data_f_grid.) Return also in this case.
  if (f_extent < 1){
  }

  // This is the part of f_grid that lies inside the spectrum data band
  const ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

  // Determine the index of the first grid points in the spectrum band that 
  // lies outside or on the boundary of the part of f_grid that lies inside
  // the spectral band.
  const Numeric f_grid_fmin = f_grid[i_fstart];
  const Numeric f_grid_fmax = f_grid[i_fstop];

  const Numeric* data_f_grid_begin = data_f_grid.unsafe_data_handle();
  const Numeric* data_f_grid_end = data_f_grid_begin + data_f_grid.size() - 1;
  const Index i_data_fstart =
      std::distance(
          data_f_grid_begin,
          std::upper_bound(data_f_grid_begin, data_f_grid_end, f_grid_fmin)) -
      1;
  const Index i_data_fstop = std::distance(
      data_f_grid_begin,
      std::upper_bound(
          data_f_grid_begin + i_data_fstart, data_f_grid_end, f_grid_fmax));

  // Extent for active data frequency vector:
  const Index data_f_extent = i_data_fstop - i_data_fstart + 1;

  // For this range the interpolation will find place.
  const Range active_range(i_data_fstart, data_f_extent);
  const ConstVectorView data_f_grid_active = data_f_grid[active_range];

  // Check if frequency is inside the range covered by the data:
  chk_interpolation_grids("Frequency interpolation for cross sections",
                          data_f_grid,
                          f_grid_active);

  {
    // Find frequency grid positions:
    ArrayOfGridPos f_gp(f_grid_active.nelem());
    gridpos(f_gp, data_f_grid_active, f_grid_active, 0);

    Matrix itw(f_gp.nelem(), 2);
    interpweights(itw, f_gp);

    for(int i=0; i < 4; i++){
      interp(int_data(Range(i_fstart, f_extent),i), itw, 
        sun_spectrum_raw.data(active_range, i), f_gp);
    }
  }  
  
  // Padding
  if (f_extent < nf){
    if (temperature == -1){
      ARTS_USER_ERROR(
      "f_grid is (partially) outside the sun spectrum data"
      "Set temperature to zero to have a padding of "
      "0 or a fitting effective temperature"
      "For further information take a look at the "
      "documentation for regrid_sun_spectrum")
    }
    if (temperature > 0){
      for (int i=0; i < i_fstart; i++){
        int_data(i,0) = planck(f_grid[i], temperature);
      }
      for (Index i=f_extent; i < nf; i++){
        int_data(i,0) = planck(f_grid[i], temperature);
      }
    }
  }

  return int_data;
}


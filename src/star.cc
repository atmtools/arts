/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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
/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   star.cc
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Functions needed for radiative transfer with direct sources.
*/

#include "star.h"
#include "auto_md.h"
#include "agenda_class.h"
#include "debug.h"
#include "matpack.h"
#include "matpackI.h"
#include "messages.h"
#include "physics_funcs.h"
#include "propagationmatrix.h"
#include "geodetic.h"
#include "arts.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric sigma = 5.670374419184429453970e-8; // Stefan–Boltzmann constant W m-2 K-4

/*===========================================================================
  === The functions
  ===========================================================================*/

std::ostream& operator<<(std::ostream& os, const Star& star) {
  os << "Star: " << star.description;
  os << " Radius: " << star.radius << "m ";
  os << " Distance: " << star.distance << "m \n";
  os << " Latitude: " << star.latitude << "° \n";
  os << " Longitude: " << star.longitude << "° \n";
  os << " Spectrum [W/m2/Hz]: \n" << star.spectrum ;
  return os;
}

void get_scattered_starsource(Workspace& ws,
                              RadiationVector& scattered_starlight,
                              ArrayOfRadiationVector& dscattered_starlight,
                              const Vector& f_grid,
                              const Numeric& p,
                              const Numeric& T,
                              const Vector& vmr,
                              const Matrix& transmitted_starlight,
                              const Vector& in_los,
                              const Vector& out_los,
                              const Agenda& gas_scattering_agenda) {
  PropagationMatrix K_sca;
  TransmissionMatrix sca_mat;
  Vector sca_fct_dummy;

  // calculate gas scattering properties
  gas_scattering_agendaExecute(ws,
                               K_sca,
                               sca_mat,
                               sca_fct_dummy,
                               f_grid,
                               p,
                               T,
                               vmr,
                               in_los,
                               out_los,
                               0,
                               gas_scattering_agenda);

  //some basic quantities
  Index ns = transmitted_starlight.ncols();
  Index nf = f_grid.nelem();

  //allocate and resize
  RadiationVector scattered_starlight_temp(1, ns);

  Matrix mat_temp(1, ns,0.);
  // Calculate the scattered radiation
  for (Index i_f = 0; i_f < nf; i_f++) {
    mat_temp(0,joker) =  transmitted_starlight(i_f, joker);
    scattered_starlight_temp = mat_temp;//transmitted_starlight(i_f, joker);
    scattered_starlight_temp.leftMul(sca_mat);

    for (Index j = 0; j < ns; j++) {
      scattered_starlight(i_f, j) =
          scattered_starlight_temp(0, j) * K_sca.Kjj(0, 0)[i_f];
    }
  }
}

void get_star_background(Matrix& iy,
                         Index& stars_visible,
                         const ArrayOfStar& stars,
                         const Ppath& ppath,
                         const Index& atmosphere_dim,
                         const Vector& refellipsoid) {
  const Index np = ppath.np;

  //set visibilty flag to default
  stars_visible = 0;

  Vector rtp_pos, rtp_los;
  rtp_pos.resize(atmosphere_dim);
  rtp_pos = ppath.pos(np - 1, Range(0, atmosphere_dim));
  rtp_los.resize(ppath.los.ncols());
  rtp_los = ppath.los(np - 1, joker);

  for (Index i_star = 0; i_star < stars.nelem(); i_star++) {
    get_star_radiation(iy, stars_visible, stars[i_star], rtp_pos, rtp_los, refellipsoid);
  }
}

void get_star_radiation(Matrix& iy,
                        Index& stars_visible,
                         const Star& stars,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& refellipsoid) {

  //Calculate earth centric radial component of star_pos and rtp_pos.
  const Numeric R_star = stars.distance;
  const Numeric R_rte = refell2r(refellipsoid, rtp_pos[1]) + rtp_pos[0];

  //Transform to cartesian coordinate system
  Numeric r_star_x, r_star_y, r_star_z;
  Numeric r_rte_x, r_rte_y, r_rte_z;
  Numeric r_los_x, r_los_y, r_los_z;

  // r_star
  sph2cart(r_star_x, r_star_y, r_star_z, R_star, stars.latitude, stars.longitude);

  // r_rte, r_los
  poslos2cart(r_rte_x,
              r_rte_y,
              r_rte_z,
              r_los_x,
              r_los_y,
              r_los_z,
              R_rte,
              rtp_pos[1],
              rtp_pos[2],
              rtp_los[0],
              rtp_los[1]);

  //Calculate vector of line of sight and unit vector pointing from
  //ppath point to the star.
  const Numeric r_ps_x = r_star_x - r_rte_x;
  const Numeric r_ps_y = r_star_y - r_rte_y;
  const Numeric r_ps_z = r_star_z - r_rte_z;

  //abs value of r_ps
  const Numeric r_ps =
      sqrt(r_ps_x * r_ps_x + r_ps_y * r_ps_y + r_ps_z * r_ps_z);

  //abs value of r_los
  const Numeric r_glos =
      sqrt(r_los_x * r_los_x + r_los_y * r_los_y + r_los_z * r_los_z);

  //Calculate angle beta between line of sight and the line between ppath point and the star
  //using scalar product
  const Numeric cos_beta =
      (r_ps_x * r_los_x + r_ps_y * r_los_y + r_ps_z * r_los_z) / (r_ps * r_glos);
  const Numeric beta = acos(cos_beta);

  // angular radius of star
  const Numeric alpha = atan(stars.radius / r_ps);

  //Check if we see the star. We see the star if the angle beta is smaller than
  // the angular radius alpha of the star.
  if (beta <= alpha) {
    //Here we assume that the star radiates isotropically.
    Matrix star_radiance = stars.spectrum;
    star_radiance /= PI;

    iy += star_radiance;

    // visibility flag
    stars_visible = 1;
  }
}

void get_transmitted_starlight(
    Workspace& ws,
    Matrix& transmitted_starlight,
    ArrayOfTensor3 dtransmitted_starlight,
    Vector& star_rte_los,
    Index& star_path_ok,
    const Vector& rte_pos,
    const Index& i_star,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
    const Tensor3& t_field,
    const EnergyLevelMap& nlte_field,
    const Tensor4& vmr_field,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Tensor3& wind_u_field,
    const Tensor3& wind_v_field,
    const Tensor3& wind_w_field,
    const Tensor3& mag_u_field,
    const Tensor3& mag_v_field,
    const Tensor3& mag_w_field,
    const Matrix& z_surface,
    const Vector& refellipsoid,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Index& gas_scattering_do,
    const Index& jacobian_do,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfStar& stars,
    const Numeric& rte_alonglos_v,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& water_p_eq_agenda,
    const Agenda& gas_scattering_agenda,
    const Agenda& ppath_step_agenda,
    const Verbosity& verbosity) {
  Vector star_pos(3);
  Ppath star_ppath;
  Numeric ppath_lraytrace2 = ppath_lraytrace;

  //dummy variables needed for the output and input of iyTransmission and
  // gas_scattering_agenda
  ArrayOfMatrix iy_aux_dummy;
  ArrayOfString iy_aux_vars_dummy;
  Vector ppvar_p_dummy;
  Vector ppvar_t_dummy;
  EnergyLevelMap ppvar_nlte_dummy;
  Matrix ppvar_vmr_dummy;
  Matrix ppvar_wind_dummy;
  Matrix ppvar_mag_dummy;
  Matrix ppvar_pnd_dummy;
  Matrix ppvar_f_dummy;
  Tensor3 ppvar_iy_dummy;
  Tensor4 ppvar_trans_cumulat_dummy;
  Tensor4 ppvar_trans_partial_dummy;
  const Tensor4 pnd_field_dummy;
  const ArrayOfTensor4 dpnd_field_dx_dummy;
  const ArrayOfString scat_species_dummy;
  const ArrayOfArrayOfSingleScatteringData scat_data_dummy;
  TransmissionMatrix sca_mat_dummy;
  Vector sca_fct_dummy;
  const Vector in_los_dummy;
  const Vector out_los_dummy;

  star_pos = {
      stars[i_star].distance, stars[i_star].latitude, stars[i_star].longitude};

  // we need the distance to the star relative to the surface
  star_pos[0] =
      star_pos[0] -
      pos2refell_r(atmosphere_dim, refellipsoid, lat_grid, lon_grid, star_pos);

  // get the line of sight direction from star i_star to ppath point
  rte_losGeometricFromRtePosToRtePos2(star_rte_los,
                                      atmosphere_dim,
                                      lat_grid,
                                      lon_grid,
                                      refellipsoid,
                                      rte_pos,
                                      star_pos,
                                      verbosity);

  // calculate ppath (star path) from star to ppath point
  ppathFromRtePos2(ws,
                   star_ppath,
                   star_rte_los,
                   ppath_lraytrace2,
                   ppath_step_agenda,
                   atmosphere_dim,
                   p_grid,
                   lat_grid,
                   lon_grid,
                   z_field,
                   f_grid,
                   refellipsoid,
                   z_surface,
                   rte_pos,
                   star_pos,
                   ppath_lmax,
                   2e-5,
                   5.,
                   0.5,
                   verbosity);

  // Check that star path hast more than 1 ppath points and that space is background
  star_path_ok = star_ppath.np > 1 && ppath_what_background(star_ppath) == 9;

  if (star_path_ok) {
    //get the TOA distance to earth center.
    Numeric R_TOA =
        refell2r(refellipsoid, star_ppath.pos(star_ppath.np - 1, 1)) +
        star_ppath.pos(star_ppath.np - 1, 0);

    //get the distance between sun and star_ppath at TOA
    Numeric R_Star2Toa;
    distance3D(R_Star2Toa,
               R_TOA,
               star_ppath.pos(star_ppath.np - 1, 1),
               star_ppath.pos(star_ppath.np - 1, 2),
               star_pos[0],
               star_pos[1],
               star_pos[2]);

    // Scale the incoming star_irradiance spectrum
    Matrix star_spectrum = stars[i_star].spectrum;
    star_spectrum *= stars[i_star].radius * stars[i_star].radius;
    star_spectrum /=
        (stars[i_star].radius * stars[i_star].radius + R_Star2Toa * R_Star2Toa);

    // calculate transmission through atmosphere
    iyTransmissionStandard(ws,
                           transmitted_starlight,
                           iy_aux_dummy,
                           dtransmitted_starlight,
                           ppvar_p_dummy,
                           ppvar_t_dummy,
                           ppvar_nlte_dummy,
                           ppvar_vmr_dummy,
                           ppvar_wind_dummy,
                           ppvar_mag_dummy,
                           ppvar_pnd_dummy,
                           ppvar_f_dummy,
                           ppvar_iy_dummy,
                           ppvar_trans_cumulat_dummy,
                           ppvar_trans_partial_dummy,
                           stokes_dim,
                           f_grid,
                           atmosphere_dim,
                           p_grid,
                           t_field,
                           nlte_field,
                           vmr_field,
                           abs_species,
                           wind_u_field,
                           wind_v_field,
                           wind_w_field,
                           mag_u_field,
                           mag_v_field,
                           mag_w_field,
                           cloudbox_on,
                           cloudbox_limits,
                           gas_scattering_do,
                           pnd_field_dummy,
                           dpnd_field_dx_dummy,
                           scat_species_dummy,
                           scat_data_dummy,
                           iy_aux_vars_dummy,
                           jacobian_do,
                           jacobian_quantities,
                           star_ppath,
                           star_spectrum,
                           propmat_clearsky_agenda,
                           water_p_eq_agenda,
                           gas_scattering_agenda,
                           1,
                           Tensor3(),
                           rte_alonglos_v,
                           verbosity);
  }
  else{
    transmitted_starlight=stars[i_star].spectrum;
    transmitted_starlight*=0.;
  }
}

Matrix regrid_star_spectrum(const GriddedField2 &star_spectrum_raw,
                          const Vector &f_grid,
                          const Index &stokes_dim,
                          const Numeric &temperature,                        
                          const Verbosity &verbosity){
  CREATE_OUTS;
  
  const Index nf = f_grid.nelem();
  const Vector data_f_grid = star_spectrum_raw.get_numeric_grid(0);
  const Numeric data_fmin = data_f_grid[0];
  const Numeric data_fmax = data_f_grid[data_f_grid.nelem() - 1];

    // Result array
  Matrix int_data(f_grid.nelem(), stokes_dim, 0.);
  
  if (out3.sufficient_priority()) {
    ostringstream os;
    os << "    f_grid:      " << f_grid[0] << " - " << f_grid[nf - 1]
      << " Hz\n"
      << "    data_f_grid: " << data_fmin << " - " << data_fmax << " Hz\n";
    out3 << os.str();
  }

  const Numeric* f_grid_begin = f_grid.get_c_array();
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
    out3 << "All frequencies are below or above the star spectrum data";
  }

  const Index f_extent = i_fstop - i_fstart + 1;

  if (out3.sufficient_priority()) {
    ostringstream os;
    os << "    " << f_extent << " frequency extraction points starting at "
        << "frequency index " << i_fstart << ".\n";
    out3 << os.str();
  }

  // If f_extent is less than one, then the entire data_f_grid is between two
  // grid points of f_grid. (So that we do not have any f_grid points inside
  // data_f_grid.) Return also in this case.
  if (f_extent < 1){
    out3 << "The entire data lies inbetween two f_grid points (So no f_grid"
        << " point inside the data_f_grid)";
  }

  // This is the part of f_grid that lies inside the spectrum data band
  const ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

  // Determine the index of the first grid points in the spectrum band that 
  // lies outside or on the boundary of the part of f_grid that lies inside
  // the spectral band.
  const Numeric f_grid_fmin = f_grid[i_fstart];
  const Numeric f_grid_fmax = f_grid[i_fstop];

  const Numeric* data_f_grid_begin = data_f_grid.get_c_array();
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

    for(int i=0; i < stokes_dim; i++){
      interp(int_data(Range(i_fstart, f_extent),i), itw, 
        star_spectrum_raw.data(active_range, i), f_gp);
    }
  }  
  
  // Padding
  if (i_fstart < 0 || f_extent == nf){
    if (temperature == -1){
      ARTS_USER_ERROR(
      "f_grid is (partially) outside the star spectrum data"
      "Set temperature to zero to have a padding of "
      "0 or a fitting effective temperature"
      "For further information take a look at the "
      "documentation for regrid_star_spectrum")
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


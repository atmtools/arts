/**
 * @file   mc_antenna.cc
 * @author Cory Davis <cdavis@staffmail.ed.ac.uk>
 * @date   2005-12-01
 *
 * @brief  Monte Carlo Antenna implementation
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "mc_antenna.h"

#include <arts_conversions.h>

#include <cfloat>

Matrix33 rotmat_enu(const Vector2& prop_los) {
  using Conversion::cosd, Conversion::sind;

  const Numeric cza = cosd(prop_los[0]);
  const Numeric sza = sind(prop_los[0]);
  const Numeric caa = cosd(prop_los[1]);
  const Numeric saa = sind(prop_los[1]);

  Matrix33 R_ant2enu{};

  R_ant2enu[0, 0] = -cza * saa;
  R_ant2enu[0, 1] = caa;
  R_ant2enu[0, 2] = sza * saa;

  R_ant2enu[1, 0] = -cza * caa;
  R_ant2enu[1, 1] = -saa;
  R_ant2enu[1, 2] = sza * caa;

  R_ant2enu[2, 0] = sza;
  R_ant2enu[2, 1] = 0.0;
  R_ant2enu[2, 2] = cza;

  return R_ant2enu;
}

Muelmat rotmat_stokes(Numeric f1_dir,
                      Numeric f2_dir,
                      const Matrix33& R_f1,
                      const Matrix33& R_f2) {
  const Numeric flip = f1_dir * f2_dir;
  Numeric cos_pra1, sin_pra1, cos_pra2, sin_pra2;

  cos_pra1 = dot(R_f1[joker, 0], R_f2[joker, 0]);
  sin_pra1 = f2_dir * dot(R_f1[joker, 0], R_f2[joker, 1]);
  sin_pra2 = f1_dir * dot(R_f1[joker, 1], R_f2[joker, 0]);
  cos_pra2 = f1_dir * f2_dir * dot(R_f1[joker, 1], R_f2[joker, 1]);

  Muelmat R_pra{};

  R_pra[0, 0] = 1.0;
  R_pra[1, 1] = 2 * cos_pra1 * cos_pra1 - 1.0;
  R_pra[1, 2] = flip * 2 * cos_pra1 * sin_pra1;
  R_pra[2, 1] = 2 * cos_pra2 * sin_pra2;
  R_pra[2, 2] = flip * (2 * cos_pra2 * cos_pra2 - 1.0);
  R_pra[3, 3] = flip * 1.0;

  return R_pra;
}

void MCAntenna::set_pencil_beam() { atype = AntennaType::PencilBeam; }

void MCAntenna::set_gaussian(const Numeric& za_sigma, const Numeric& aa_sigma) {
  atype    = AntennaType::Gaussian;
  sigma_za = za_sigma;
  sigma_aa = aa_sigma;
}

void MCAntenna::set_gaussian_fwhm(const Numeric& za_fwhm,
                                  const Numeric& aa_fwhm) {
  atype    = AntennaType::Gaussian;
  sigma_za = za_fwhm / 2.3548;
  sigma_aa = aa_fwhm / 2.3548;
}

void MCAntenna::set_lookup(ConstVectorView za_grid_,
                           ConstVectorView aa_grid_,
                           ConstMatrixView G_lookup_) {
  atype    = AntennaType::Lookup;
  za_grid  = za_grid_;
  aa_grid  = aa_grid_;
  G_lookup = G_lookup_;
}

Numeric MCAntenna::return_los(const Matrix33& R_return,
                              const Matrix33& R_enu2ant) const {
  using Conversion::atan2d;

  Numeric z, term_el, term_az;
  Numeric ant_el, ant_az;
  Vector3 k_vhk;

  switch (atype) {
    case AntennaType::PencilBeam: return 1.0;

    case AntennaType::Gaussian:
      mult(k_vhk, R_enu2ant, R_return[joker, 2]);

      // Assume Gaussian is narrow enough that response is 0 beyond 90 degrees
      // Same assumption is made for drawing samples (draw_los)
      if (k_vhk[2] > 0) {
        ant_el  = atan2d(k_vhk[0], k_vhk[2]);
        ant_az  = atan2d(k_vhk[1], k_vhk[2]);
        term_el = ant_el / sigma_za;
        term_az = ant_az / sigma_aa;
        z       = term_el * term_el + term_az * term_az;
        return std::exp(-0.5 * z);
      }

      return 0.0;

    default: ARTS_USER_ERROR("invalid Antenna type.")
  }
}

std::pair<Vector2, Matrix33> MCAntenna::draw_los(
    RandomNumberGenerator<>& rng,
    const Matrix33& R_ant2enu,
    const Vector2& bore_sight_los) const {
  using Conversion::tand, Conversion::atan2d, Conversion::acosd;

  Numeric ant_el, ant_az, ant_r;
  Numeric tel, taz;
  Vector3 k_vhk;

  std::pair<Vector2, Matrix33> result;
  auto& [sampled_rte_los, R_los] = result;

  switch (atype) {
    case AntennaType::PencilBeam:
      sampled_rte_los = bore_sight_los;
      R_los           = R_ant2enu;
      break;

    case AntennaType::Gaussian:

      ant_el = 91;
      ant_az = 91;

      // Assume Gaussian is narrow enough that response is 0 beyond 90 degrees
      // Same assumption is made for radar return samples (return_los)
      {
        auto za_sample = rng.get<std::normal_distribution>(0.0, sigma_za);
        while (ant_el >= 90) {
          ant_el = za_sample();
        }
        auto aa_sample = rng.get<std::normal_distribution>(0.0, sigma_aa);
        while (ant_az >= 90) {
          ant_az = aa_sample();
        }
      }

      // Propagation direction
      tel      = tand(ant_el);
      taz      = tand(ant_az);
      ant_r    = std::sqrt(1 + tel * tel + taz * taz);
      k_vhk[0] = tel / ant_r;
      k_vhk[1] = taz / ant_r;
      k_vhk[2] = 1.0 / ant_r;
      mult(R_los[joker, 2], R_ant2enu, k_vhk);

      sampled_rte_los[0] = acosd(R_los[2, 2]);

      // Horizontal polarization basis
      // If drawn los is at zenith or nadir, assume same azimuth as boresight
      if ((1.0 - std::abs(R_los[2, 2])) < DBL_EPSILON) {
        // H is aligned with H of bs, use row not column because tranpose
        R_los[joker, 1]    = R_ant2enu[1, joker];
        sampled_rte_los[1] = bore_sight_los[1];
      } else {
        const Vector uhat{0.0, 0.0, 1.0};
        Numeric magh;
        sampled_rte_los[1] = atan2d(R_los[0, 2], R_los[1, 2]);
        cross3(R_los[joker, 1], R_los[joker, 2], uhat);
        magh             = hypot(R_los[joker, 1]);
        R_los[joker, 1] /= magh;
      }

      // Vertical polarization basis
      cross3(R_los[joker, 0], R_los[joker, 1], R_los[joker, 2]);

      break;

    default: ARTS_USER_ERROR("invalid Antenna type.")
  }

  return result;
}

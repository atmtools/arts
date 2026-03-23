/**
  @file   m_radar.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @author Richard Larsson <richard.larsson@ri.se>
  @date   2025-03-12

  @brief  Workspace functions related to simulation of cloud radars.

  Ported from ARTS2 m_cloudradar.cc to ARTS3 patterns.
 */

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <debug.h>
#include <mie.h>
#include <path_point.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <complex>

namespace {

/** Compute conversion factors from backscatter coefficient to Ze.
 *
 * fac = (4e18 / pi^4) * lambda^4 / |K|^2
 *
 * If k2 <= 0, the Liebe93 model computes |K|^2 at the reference temperature.
 */
Vector ze_cfac_calc(const AscendingGrid& f_grid,
                    const Numeric ze_tref,
                    const Numeric k2) {
  const Index nf = f_grid.size();
  Vector fac(nf);

  constexpr Numeric a = 4e18 / Math::pow4(Constant::pi);

  for (Index iv = 0; iv < nf; iv++) {
    Numeric K2;
    if (k2 > 0) {
      K2 = k2;
    } else {
      auto n = scattering::refr_index_water_liebe93(f_grid[iv], ze_tref);
      K2     = std::norm(scattering::dielectric_factor(n));
    }

    const Numeric la = Constant::speed_of_light / f_grid[iv];
    fac[iv]          = a * Math::pow4(la) / K2;
  }

  return fac;
}

}  // namespace

/* Workspace method: Doxygen documentation will be auto-generated */
void radar_spectral_radSingleScat(
    StokvecMatrix& radar_spectral_rad,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& freq_grid,
    const SurfaceField& surf_field,
    const ArrayOfPropmatVector& spectral_propmat_path,
    const ArrayOfMuelmatVector& radar_bulk_backscatter,
    const StokvecVector& radar_spectral_rad_transmitter,
    const Numeric& radar_pext_scaling) try {
  ARTS_TIME_REPORT

  const Size nf = freq_grid.size();
  const Size np = ray_path.size();

  ARTS_USER_ERROR_IF(np < 2, "Need at least two path points, got {}", np);
  ARTS_USER_ERROR_IF(spectral_propmat_path.size() != np,
                     "spectral_propmat_path size ({}) != ray_path size ({})",
                     spectral_propmat_path.size(),
                     np);
  ARTS_USER_ERROR_IF(
      radar_spectral_rad_transmitter.size() != nf,
      "radar_spectral_rad_transmitter size ({}) != freq_grid size ({})",
      radar_spectral_rad_transmitter.size(),
      nf);
  ARTS_USER_ERROR_IF(radar_pext_scaling < 0 || radar_pext_scaling > 2,
                     "radar_pext_scaling must be in [0, 2], got {}",
                     radar_pext_scaling);
  ARTS_USER_ERROR_IF(radar_bulk_backscatter.size() != np,
                     "radar_bulk_backscatter size ({}) != ray_path size ({})",
                     radar_bulk_backscatter.size(),
                     np);

  const Vector r = path::distance(ray_path, surf_field.ellipsoid);

  // Build per-layer transmission matrices.  T[0] is identity.
  using rtepack::muelmat;
  Array<MuelmatVector> lyr_tra(np, MuelmatVector(nf, muelmat::id()));

  for (Size ip = 1; ip < np; ip++) {
    for (Size iv = 0; iv < nf; iv++) {
      const auto kavg = rtepack::avg(spectral_propmat_path[ip - 1][iv],
                                     spectral_propmat_path[ip][iv]);
      lyr_tra[ip][iv] = rtepack::exp(kavg, r[ip - 1]);
    }
  }

  // Cumulative forward transmission: PiTf[ip] = T[1]*T[2]*...*T[ip]
  // This is the one-way transmission from the sensor to path point ip.
  // For the radar return, the backscattered signal travels the same path
  // back (commutative/reciprocal case), so both legs use PiTf.
  auto PiTf = rtepack::forward_cumulative_transmission(lyr_tra);

  // Use the pre-computed bulk backscatter matrix
  const auto& Z = radar_bulk_backscatter;

  // Compute radar return: PiTf * Z * PiTf * I0
  // Two-way transmission: down from sensor to scatterer, backscatter,
  // then back up from scatterer to sensor.
  radar_spectral_rad.resize(nf, np);
  radar_spectral_rad = rtepack::stokvec{};

  for (Size ip = 0; ip < np; ip++) {
    for (Size iv = 0; iv < nf; iv++) {
      radar_spectral_rad[iv, ip] = PiTf[ip][iv] * Z[ip][iv] * PiTf[ip][iv] *
                                   radar_spectral_rad_transmitter[iv];
    }
  }
}
ARTS_METHOD_ERROR_CATCH

/* Workspace method: Doxygen documentation will be auto-generated */
void measurement_vecFromRadarSpectralRad(
    Vector& measurement_vec,
    const StokvecMatrix& radar_spectral_rad,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& freq_grid,
    const SurfaceField& surf_field,
    const Vector& radar_range_bins,
    const String& radar_unit,
    const Numeric& radar_ze_tref,
    const Numeric& radar_k2,
    const Numeric& radar_dbze_min) try {
  ARTS_TIME_REPORT

  const Index nf   = freq_grid.size();
  const Index np   = static_cast<Index>(ray_path.size());
  const Index nbin = radar_range_bins.size() - 1;

  ARTS_USER_ERROR_IF(np < 2, "Need at least two path points");
  ARTS_USER_ERROR_IF(
      radar_spectral_rad.nrows() != nf || radar_spectral_rad.ncols() != np,
      "radar_spectral_rad must have shape [nf, np], got [{}, {}]",
      radar_spectral_rad.nrows(),
      radar_spectral_rad.ncols());
  ARTS_USER_ERROR_IF(nbin < 1, "Need at least two range bin edges");
  ARTS_USER_ERROR_IF(!is_increasing(radar_range_bins),
                     "radar_range_bins must be strictly increasing");

  // Determine if range bins are altitude [m] (large values) or time [s]
  const bool is_z = max(radar_range_bins) > 1.0;

  ARTS_USER_ERROR_IF(!is_z && min(radar_range_bins) < 0,
                     "Time-based range bins must be non-negative");

  // Conversion factors for Ze
  Vector cfac(nf, 1.0);
  Numeric ze_min = 0;

  const bool unit_is_ze   = (radar_unit == "Ze");
  const bool unit_is_dbze = (radar_unit == "dBZe");

  ARTS_USER_ERROR_IF(
      radar_unit != "1" && !unit_is_ze && !unit_is_dbze,
      "radar_unit must be \"1\", \"Ze\", or \"dBZe\", got \"{}\"",
      radar_unit);

  if (unit_is_ze || unit_is_dbze) {
    cfac = ze_cfac_calc(freq_grid, radar_ze_tref, radar_k2);
  }
  if (unit_is_dbze) {
    ze_min = std::pow(10.0, radar_dbze_min / 10.0);
  }

  // Compute range along path (altitude or round-trip time)
  Vector range(np);
  if (is_z) {
    for (Index ip = 0; ip < np; ip++) {
      range[ip] = ray_path[ip].altitude();
    }
  } else {
    const Vector r = path::distance(ray_path, surf_field.ellipsoid);
    range[0]       = 0;
    for (Index ip = 1; ip < np; ip++) {
      const Numeric ngroup_avg =
          0.5 * (ray_path[ip - 1].ngroup + ray_path[ip].ngroup);
      range[ip] = range[ip - 1] +
                  2.0 * r[ip - 1] * ngroup_avg / Constant::speed_of_light;
    }
  }

  const Numeric range_end1 = std::min(range[0], range[np - 1]);
  const Numeric range_end2 = std::max(range[0], range[np - 1]);

  // Output: one value per frequency per bin
  const Index ntot = nf * nbin;
  measurement_vec.resize(ntot);
  measurement_vec = std::numeric_limits<Numeric>::quiet_NaN();

  for (Index b = 0; b < nbin; b++) {
    if (radar_range_bins[b] >= range_end2 ||
        radar_range_bins[b + 1] <= range_end1)
      continue;

    const Numeric blim1  = std::max(radar_range_bins[b], range_end1);
    const Numeric blim2  = std::min(radar_range_bins[b + 1], range_end2);
    const Numeric bwidth = blim2 - blim1;
    if (bwidth <= 0) continue;

    // Trapezoidal weights for averaging over the bin
    Vector hbin(np, 0.0);
    for (Index ip = 0; ip < np; ip++) {
      Numeric seg_lo, seg_hi;
      if (ip == 0) {
        seg_lo = range[0];
        seg_hi = 0.5 * (range[0] + range[1]);
      } else if (ip == np - 1) {
        seg_lo = 0.5 * (range[np - 2] + range[np - 1]);
        seg_hi = range[np - 1];
      } else {
        seg_lo = 0.5 * (range[ip - 1] + range[ip]);
        seg_hi = 0.5 * (range[ip] + range[ip + 1]);
      }
      if (seg_lo > seg_hi) std::swap(seg_lo, seg_hi);

      const Numeric olo = std::max(seg_lo, blim1);
      const Numeric ohi = std::min(seg_hi, blim2);
      if (ohi > olo) hbin[ip] = (ohi - olo) / bwidth;
    }

    for (Index iv = 0; iv < nf; iv++) {
      const Index iout = iv * nbin + b;

      Numeric val = 0;
      for (Index ip = 0; ip < np; ip++) {
        val += hbin[ip] * radar_spectral_rad[iv, ip].I();
      }
      val *= cfac[iv];

      if (unit_is_dbze) {
        measurement_vec[iout] =
            val <= ze_min ? radar_dbze_min : 10.0 * std::log10(val);
      } else {
        measurement_vec[iout] = val;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

/* Workspace method: Doxygen documentation will be auto-generated */
void radar_ze_cfacCalc(Vector& radar_ze_cfac,
                       const AscendingGrid& freq_grid,
                       const Numeric& radar_ze_tref,
                       const Numeric& radar_k2) try {
  ARTS_TIME_REPORT

  radar_ze_cfac = ze_cfac_calc(freq_grid, radar_ze_tref, radar_k2);
}
ARTS_METHOD_ERROR_CATCH

/* Workspace method: Doxygen documentation will be auto-generated */
void radar_bulk_backscatterFromScat(
    ArrayOfMuelmatVector& radar_bulk_backscatter,
    ArrayOfPropmatVector& spectral_propmat_path,
    const ArrayOfScatteringSpecies& scat_species,
    const ArrayOfAtmPoint& atm_path,
    const AscendingGrid& freq_grid,
    const Numeric& radar_pext_scaling) try {
  ARTS_TIME_REPORT

  const Index nf = freq_grid.size();
  const Index np = static_cast<Index>(atm_path.size());

  ARTS_USER_ERROR_IF(np < 2, "Need at least two path points, got {}", np);

  // Determine if we already have gas absorption in spectral_propmat_path
  const bool have_existing_propmat =
      static_cast<Index>(spectral_propmat_path.size()) == np;

  if (!have_existing_propmat) {
    spectral_propmat_path.resize(np);
    for (Index ip = 0; ip < np; ip++) {
      spectral_propmat_path[ip] = PropmatVector(nf, rtepack::propmat{});
    }
  }

  radar_bulk_backscatter.resize(np);
  for (Index ip = 0; ip < np; ip++) {
    radar_bulk_backscatter[ip] = MuelmatVector(nf, rtepack::muelmat{0.0});
  }

  // For each path point, accumulate radar backscatter and extinction
  // from all scattering species via the variant dispatch.
  for (Index ip = 0; ip < np; ip++) {
    const auto& atm = atm_path[ip];

    for (const auto& spec : scat_species.species) {
      MuelmatVector Z_back;
      PropmatVector K_ext;

      std::visit(
          [&](const auto& s) {
            if constexpr (requires {
                            s.get_radar_single_scat(
                                Z_back, K_ext, atm, freq_grid);
                          }) {
              s.get_radar_single_scat(Z_back, K_ext, atm, freq_grid);

              for (Index iv = 0; iv < nf; iv++) {
                radar_bulk_backscatter[ip][iv] += Z_back[iv];
                spectral_propmat_path[ip][iv].A() +=
                    K_ext[iv].A() * radar_pext_scaling;
              }
            }
          },
          spec);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

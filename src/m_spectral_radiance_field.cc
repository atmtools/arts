#include <arts_omp.h>
#include <workspace.h>

#include "math_funcs.h"

namespace {
auto ray_path_propagation_matrixProfile(
    const Workspace& ws,
    const Agenda& propagation_matrix_agenda,
    const AscendingGrid& frequency_grid,
    const ArrayOfAtmPoint& ray_path_atmospheric_point) {
  struct ray_path_propagation_matrixFromPathOutput {
    ArrayOfPropmatVector k;
    ArrayOfStokvecVector s;
    ArrayOfPropmatMatrix dk;
    ArrayOfStokvecMatrix ds;
  };
  ray_path_propagation_matrixFromPathOutput out;

  const Size np = ray_path_atmospheric_point.size();
  out.k.resize(np);
  out.s.resize(np);
  out.dk.resize(np);
  out.ds.resize(np);

  String error{};

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    try {
      propagation_matrix_agendaExecute(ws,
                                       out.k[ip],
                                       out.s[ip],
                                       out.dk[ip],
                                       out.ds[ip],
                                       frequency_grid,
                                       {},
                                       {},
                                       {},
                                       {},
                                       ray_path_atmospheric_point[ip],
                                       propagation_matrix_agenda);
    } catch (const std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  return out;
}

StokvecVector interp(const StokvecConstMatrixView& data,
                     const Vector& zenith_grid,
                     const Numeric& zenith) {
  const Size NF = data.shape().back();

  const FixedLagrangeInterpolation<1> za(0, zenith, zenith_grid);

  StokvecVector out(NF);

  for (Size i = 0; i < NF; i++) out[i] = interp(data[joker, i], za);

  return out;
}

const Vector zenith_level_limb(const AscendingGrid& altitude_grid,
                               const Vector2& ellipsoid,
                               const Numeric& latitude,
                               const Numeric& longitude,
                               const Numeric& azimuth) {
  const Size NA = altitude_grid.size();
  Vector zenith_limb(NA, 90.0);

  for (Size i = 1; i < NA; i++) {
    const Numeric alt_low  = altitude_grid[i - 1];
    const Numeric alt_this = altitude_grid[i];
    const Vector3 pos{alt_this, latitude, longitude};
    zenith_limb[i] =
        path::geometric_tangent_zenith(pos, ellipsoid, alt_low, azimuth);
  }

  return zenith_limb;
}
}  // namespace

void zenith_gridProfilePseudo2D(AscendingGrid& zenith_grid,
                                const SurfaceField& surface_field,
                                const AscendingGrid& altitude_grid,
                                const Numeric& latitude,
                                const Numeric& longitude,
                                const Numeric& dza,
                                const Numeric& azimuth,
                                const Index& consider_limb) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(dza <= 0.0 or dza >= 180.0,
                     "Delta zenith angle must be (0, 180). Given dza: {}",
                     dza);

  Vector out;

  // Helper to push into the output
  const auto append = [&](auto&& x) {
    for (auto&& za : x) out.push_back(za);
  };

  if (static_cast<bool>(consider_limb)) {
    Vector zenith_limb = zenith_level_limb(
        altitude_grid, surface_field.ellipsoid, latitude, longitude, azimuth);
    stdr::sort(zenith_limb);
    const auto [it, _] = stdr::unique(zenith_limb);

    append(binary_grid(0, std::nextafter(zenith_limb.front(), 0.0), dza));
    out.push_back(90.0);
    out.push_back(std::nextafter(90., 180.0));

    Numeric past_za = zenith_limb.front();
    for (Numeric za :
         stdr::subrange(stdr::begin(zenith_limb), it) | stdv::drop(1)) {
      const auto hit = std::nextafter(za, 0.0);
      const auto mis = std::nextafter(za, 180.0);
      append(binary_grid(past_za, hit, dza) | stdv::drop(1));
      out.push_back(mis);
      past_za = za;
    }
    stdr::sort(out);

    append(binary_grid(out.back(), 180.0, dza) | stdv::drop(1));
  } else {
    out = binary_grid(0, 180, dza);
  }

  auto [it, _] = stdr::unique(out);
  auto s       = out | stdv::take(stdr::distance(out.begin(), it));
  zenith_grid  = AscendingGrid(stdr::begin(s), stdr::end(s), std::identity{});
}

void spectral_radiance_fieldProfilePseudo2D(
    const Workspace& ws,
    StokvecSortedGriddedField6& spectral_radiance_field,
    const Agenda& propagation_matrix_agenda,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const SurfaceField& surface_field,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& zenith_grid,
    const AscendingGrid& altitude_grid,
    const Numeric& latitude,
    const Numeric& longitude,
    const Numeric& azimuth) try {
  ARTS_TIME_REPORT

  constexpr Numeric minimal_r = 0.001;

  ARTS_USER_ERROR_IF(
      not arr::same_size(altitude_grid, ray_path_atmospheric_point),
      R"(Altitude grid and atmospheric point grid must have the same size

Altitude grid size:          {}
Atmospheric point grid size: {}
)",
      altitude_grid.size(),
      ray_path_atmospheric_point.size())

  ARTS_USER_ERROR_IF(zenith_grid.empty(), "Need some zenith angles")

  const Size NA = altitude_grid.size();
  const Size NZ = zenith_grid.size();
  const Size NF = frequency_grid.size();

  spectral_radiance_field.data_name = "spectral_radiance_fieldProfilePseudo2D";

  spectral_radiance_field.resize(NA, 1, 1, NZ, 1, NF);
  spectral_radiance_field.data = 0.0;

  spectral_radiance_field.grid<0>() = altitude_grid;
  spectral_radiance_field.grid<1>() = Vector{latitude};
  spectral_radiance_field.grid<2>() = Vector{longitude};
  spectral_radiance_field.grid<3>() = zenith_grid;
  spectral_radiance_field.grid<4>() = Vector{azimuth};
  spectral_radiance_field.grid<5>() = frequency_grid;

  spectral_radiance_field.gridname<0>() = "altitude";
  spectral_radiance_field.gridname<1>() = "latitude";
  spectral_radiance_field.gridname<2>() = "longitude";
  spectral_radiance_field.gridname<3>() = "zenith";
  spectral_radiance_field.gridname<4>() = "azimuth";
  spectral_radiance_field.gridname<5>() = "frequency";

  if (NA == 0) return;

  const auto propmat_data =
      ray_path_propagation_matrixProfile(ws,
                                         propagation_matrix_agenda,
                                         frequency_grid,
                                         ray_path_atmospheric_point);

  constexpr Numeric t_spac = Constant::cosmic_microwave_background_temperature;
  const Numeric t_surf = surface_field[SurfaceKey::t].at(latitude, longitude);
  const Vector2 ell    = surface_field.ellipsoid;

  const Vector zenith_limb =
      zenith_level_limb(altitude_grid, ell, latitude, longitude, azimuth);

  ARTS_USER_ERROR_IF(
      zenith_grid.front() != 0.0 or zenith_grid.back() != 180.0 or
          not stdr::contains(zenith_grid, 90.0),
      "Zenith grid must contain 0, 90, 180, beyond this it is free-form")

  auto srad = spectral_radiance_field.data.view_as(NA, NZ, NF);

  // Background radiation
  for (Size iv = 0; iv < NF; iv++) {
    srad[0, joker, iv]      = planck(frequency_grid[iv], t_surf);
    srad[NA - 1, joker, iv] = planck(frequency_grid[iv], t_spac);
  }

  const auto update = [&](Size beg, Size end, Size iz) {
    const Numeric alt_beg = altitude_grid[beg];
    const Numeric alt_end = altitude_grid[end];

    const Vector3 pos_end{alt_end, latitude, longitude};
    const Vector2 los_end{zenith_grid[iz], azimuth};

    const auto [ecef, decef] =
        path::geodetic_poslos2ecef(pos_end, path::mirror(los_end), ell);

    const auto [r0, r1] =
        path::line_ellipsoid_altitude_intersect(alt_beg, ecef, decef, ell);

    // Fix for limb, where we don't want to hit self
    const Numeric r = (r0 < minimal_r) ? r1 : r0;

    // Again, fix for limb because it might be missed for std::nextafter(90., 0);
    const Numeric za = std::isnan(r)
                           ? 90.0
                           : path::mirror(path::ecef2geodetic_poslos(
                                              ecef + r * decef, decef, ell)
                                              .second)[0];

    StokvecVector I = interp(srad[beg], zenith_grid, za);

    if (r >= minimal_r) {
      const auto& K0 = propmat_data.k[beg];
      const auto& K1 = propmat_data.k[end];
      const auto& J0 = propmat_data.s[beg];
      const auto& J1 = propmat_data.s[end];
      const auto& T0 = ray_path_atmospheric_point[beg].temperature;
      const auto& T1 = ray_path_atmospheric_point[end].temperature;

      rtepack::nlte_step(I, frequency_grid, K0, K1, J0, J1, T0, T1, r);
    }

    return I;
  };

  // Downgoing radiation
  for (Size i = NA - 1; i > 0; i--) {
#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size iz = 0; iz < NZ; iz++) {
      if (zenith_grid[iz] < 90.0) continue;
      srad[i - 1, iz] = update(i, i - 1, iz);
    }
  }

  // Upgoing radiation
  for (Size i = 0; i < NA - 1; i++) {
#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size iz = 0; iz < NZ; iz++) {
      if (zenith_grid[iz] >= 90.0) continue;
      srad[i + 1, iz] = update(i, i + 1, iz);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

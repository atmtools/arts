#include <arts_omp.h>
#include <workspace.h>

#include <stdexcept>

void atmospheric_profileExtract(ArrayOfAtmPoint &atmospheric_profile,
                                const AtmField &atmospheric_field,
                                const AscendingGrid &altitude_grid,
                                const Numeric &latitude,
                                const Numeric &longitude) {
  const Size N = altitude_grid.size();

  atmospheric_profile.clear();
  atmospheric_profile.resize(N);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < N; i++) {
    try {
      atmospheric_profile[i] =
          atmospheric_field.at(altitude_grid[i], latitude, longitude);
    } catch (const std::exception &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}

void atmospheric_profileFromGrid(ArrayOfAtmPoint &atmospheric_profile,
                                 AscendingGrid &altitude_grid,
                                 Numeric &latitude,
                                 Numeric &longitude,
                                 const AtmField &atmospheric_field,
                                 const AtmKey &key) {
  const auto &gf3 = atmospheric_field[key].get<GeodeticField3>();
  altitude_grid   = gf3.grid<0>();

  const Size N = gf3.data.size();
  ARTS_USER_ERROR_IF(not gf3.ok() or N == 0, "Not a good field for key {}", key)
  ARTS_USER_ERROR_IF(N != altitude_grid.size(), "Not a profile key {}", key)

  latitude  = gf3.grid<1>()[0];
  longitude = gf3.grid<2>()[0];

  atmospheric_profileExtract(atmospheric_profile,
                             atmospheric_field,
                             altitude_grid,
                             latitude,
                             longitude);
}

void ray_path_atmospheric_pointFromProfile(
    ArrayOfAtmPoint &ray_path_atmospheric_point,
    const ArrayOfAtmPoint &atmospheric_profile) {
  ray_path_atmospheric_point = atmospheric_profile;
}

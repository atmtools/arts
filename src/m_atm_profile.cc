#include <arts_omp.h>
#include <workspace.h>

#include <stdexcept>

void atm_profileExtract(ArrayOfAtmPoint &atm_profile,
                        const AtmField &atm_field,
                        const AscendingGrid &alt_grid,
                        const Numeric &latitude,
                        const Numeric &longitude) {
  const Size N = alt_grid.size();

  atm_profile.clear();
  atm_profile.resize(N);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < N; i++) {
    try {
      atm_profile[i] = atm_field.at(alt_grid[i], latitude, longitude);
    } catch (const std::exception &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}

void atm_profileFromGrid(ArrayOfAtmPoint &atm_profile,
                         AscendingGrid &alt_grid,
                         Numeric &latitude,
                         Numeric &longitude,
                         const AtmField &atm_field,
                         const AtmKey &key) {
  const auto &gf3 = atm_field[key].get<GeodeticField3>();
  alt_grid        = gf3.grid<0>();

  const Size N = gf3.data.size();
  ARTS_USER_ERROR_IF(not gf3.ok() or N == 0, "Not a good field for key {}", key)
  ARTS_USER_ERROR_IF(N != alt_grid.size(), "Not a profile key {}", key)

  latitude  = gf3.grid<1>()[0];
  longitude = gf3.grid<2>()[0];

  atm_profileExtract(atm_profile, atm_field, alt_grid, latitude, longitude);
}

void ray_path_atm_pointFromProfile(ArrayOfAtmPoint &ray_path_atm_point,
                                   const ArrayOfAtmPoint &atm_profile) {
  ray_path_atm_point = atm_profile;
}

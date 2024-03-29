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

#include "check_input.h"
#include "debug.h"
#include "interpolation.h"
#include "matpack_data.h"
#include "physics_funcs.h"

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

Matrix regrid_sun_spectrum(const GriddedField2 &sun_spectrum_raw,
                          const Vector &f_grid,
                          const Numeric &temperature){
  const Index nf = f_grid.size();
  const Vector& data_f_grid = sun_spectrum_raw.grid<0>();
  const Numeric data_fmin = data_f_grid[0];
  const Numeric data_fmax = data_f_grid[data_f_grid.size() - 1];

    // Result array
  Matrix int_data(f_grid.size(), 4, 0.);

  const Numeric* f_grid_begin = f_grid.unsafe_data_handle();
  const Numeric* f_grid_end = f_grid_begin + f_grid.size();
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
    ArrayOfGridPos f_gp(f_grid_active.size());
    gridpos(f_gp, data_f_grid_active, f_grid_active, 0);

    Matrix itw(f_gp.size(), 2);
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

std::ostream& operator<<(std::ostream& os, const ArrayOfSun& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

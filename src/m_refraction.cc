/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
  \file   m_refraction.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2003-01-09

  \brief  Workspace methods releated to refraction.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts_constants.h"
#include "species_tags.h"
#include "absorption.h"
#include "check_input.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "physics_funcs.h"
#include "refraction.h"
#include "special_interp.h"

inline constexpr Numeric ELECTRON_CHARGE=-Constant::elementary_charge;
inline constexpr Numeric ELECTRON_MASS=Constant::electron_mass;
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric VACUUM_PERMITTIVITY=Constant::vacuum_permittivity;
inline constexpr Numeric TORR2PA=Conversion::torr2pa(1);

/*===========================================================================
  === WSMs for refr_index_air
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void refr_index_airFreeElectrons(Numeric& refr_index_air,
                                 Numeric& refr_index_air_group,
                                 const Vector& f_grid,
                                 const ArrayOfArrayOfSpeciesTag& abs_species,
                                 const Vector& rtp_vmr,
                                 const Index& demand_vmr_value) {
  // The expression used is found in many textbooks, e.g. Rybicki and Lightman
  // (1979). Note that the refractive index corresponds to the phase velocity.

  static const Numeric k = ELECTRON_CHARGE * ELECTRON_CHARGE /
                           (VACUUM_PERMITTIVITY * ELECTRON_MASS * 4 * PI * PI);

  Numeric edensity = 0;

  Index ife = -1;
  for (Size sp = 0; sp < abs_species.size() && ife < 0; sp++) {
    if (abs_species[sp].FreeElectrons()) {
      ife = sp;
    }
  }

  if (ife < 0) {
    if (demand_vmr_value) {
      throw std::runtime_error(
          "Free electrons not found in *abs_species* and "
          "contribution to refractive index can not be calculated.");
    }
  } else {
    edensity = rtp_vmr[ife];

    if (edensity > 0) {
      // Check that lowest frequency not too low
      // Limit at 100 GHz follows Hartmann and Leitinger, Range errors due
      // to ionospheric and tropospheric effects for signal frequencies
      // above 100 HMHz, Bull. Goed., 1984.
      if (f_grid[0] < 100e6) {
        throw std::runtime_error(
            "All frequencies must be >= 100 MHz, but "
            "this is not the case.");
      }
      if (edensity * k / (f_grid[0] * f_grid[0]) > 0.25) {
        std::ostringstream os;
        os << "All frequencies must at least be twice the plasma frequency.\n"
           << "For this particular point, the plasma frequency is: "
           << sqrt(edensity * k) / 1e6 << " MHz.";
        throw std::runtime_error(os.str());
      }

      const Numeric f = (f_grid[0] + last(f_grid)) / 2.0;
      const Numeric a = edensity * k / (f * f);
      const Numeric n = sqrt(1 - a);

      refr_index_air += n - 1;
      refr_index_air_group += 1 / n - 1;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void refr_index_airInfraredEarth(Numeric& refr_index_air,
                                 Numeric& refr_index_air_group,
                                 const Numeric& rtp_pressure,
                                 const Numeric& rtp_temperature) {
  static const Numeric bn0 = 1.000272620045304;
  static const Numeric bn02 = bn0 * bn0;
  static const Numeric bk = 288.16 * (bn02 - 1.0) / (1013.25 * (bn02 + 2.0));

  // Pa -> hPa
  const Numeric n = sqrt((2.0 * bk * rtp_pressure / 100.0 + rtp_temperature) /
                         (rtp_temperature - bk * rtp_pressure / 100.0)) -
                    1;

  refr_index_air += n;
  refr_index_air_group += n;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void refr_index_airMicrowavesEarth(Numeric& refr_index_air,
                                   Numeric& refr_index_air_group,
                                   const Numeric& rtp_pressure,
                                   const Numeric& rtp_temperature,
                                   const Vector& rtp_vmr,
                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                   const Numeric& k1,
                                   const Numeric& k2,
                                   const Numeric& k3) {
  if (abs_species.size() != static_cast<Size>(rtp_vmr.size()))
    throw std::runtime_error(
        "The number of tag groups differ between "
        "*rtp_vmr* and *abs_species*.");

  Index firstH2O = find_first_species(
      abs_species, Species::fromShortName("H2O"));

  Numeric e;
  if (firstH2O < 0)
    //throw std::runtime_error(
    //   "Water vapour is a required (must be a tag group in *abs_species*)." );
    e = 0.;
  else
    e = rtp_pressure * rtp_vmr[firstH2O];

  const Numeric n =
      (k1 * (rtp_pressure - e) + (k2 + k3 / rtp_temperature) * e) /
      rtp_temperature;

  refr_index_air += n;
  refr_index_air_group += n;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void refr_index_airMicrowavesGeneral(
    Numeric& refr_index_air,
    Numeric& refr_index_air_group,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const Vector& rtp_vmr,
    const ArrayOfArrayOfSpeciesTag& abs_species) {
  //FIXME: Shall n be rescaled for sum(VMW)=1? Doing so now, but is it correct?
  //       Short sensitivity test (tropical, dry air, ~11km tanh) shows that
  //       vmr-normalized n fits significantly better (dtanh~0.5m) than
  //       non-normalized one (dtanh~5m).

  /*
   for now, hard-coding the reference refindices and refT/p. could make that
   some re-setabel parameters (like iso ratios... also regarding storing them in
   file/data struct per species)
*/
  const Numeric p0 = 760. * TORR2PA;
  const Numeric T0 = 273.15;

  // Number of refractive species:
  const Index nrs = 6;

  // This is hardwired here and quite primitive, but should do the job.
  // Set refractive index species names.
  ArrayOfString ref_spec_names(nrs);
  ref_spec_names[0] = "N2";
  ref_spec_names[1] = "O2";
  ref_spec_names[2] = "CO2";
  ref_spec_names[3] = "H2";
  ref_spec_names[4] = "He";
  ref_spec_names[5] = "H2O";

  // Set reference refractive indices
  // Values from Newell and Baird, 1965
  Vector ref_n(nrs);
  ref_n[0] = 293.81e-6;
  ref_n[1] = 266.95e-6;
  ref_n[2] = 495.16e-6;
  ref_n[3] = 135.77e-6;
  ref_n[4] = 34.51e-6;
  // value for H2O from H2O contribution according to refr_index_airMicrowavesEarth
  // at reference conditions
  // that is: n_H2O = p_ref/T_ref * (k2 + k3/Tref)
  ref_n[5] = 5338.89e-6;

  // Checks
  if (abs_species.size() != static_cast<Size>(rtp_vmr.size()))
    throw std::runtime_error(
        "The number of tag groups differ between "
        "*rtp_vmr* and *abs_species*.");
  /*
   further checks:
   ? non-neg T
   ?
*/

  // Data management
  // Find the location of all refractive species in abs_species. Set to -1 if
  // not found. The length of array ref_spec_locations is the number of
  // considered refractive species (in this method: N2, O2, CO2, H2, He).
  // The value means:
  // -1 = not in abs_species
  // N  = species is number N in abs_species

  //Can't use this one as it inside gets the broadening species names and
  //number. Also, we would have to make a workaround for this_species. So,
  //instead we use a modified version of this function directly included here.
  /*find_broad_spec_locations(ref_spec_locations,
                             abs_species,
                             this_species);*/

  ArrayOfIndex ref_spec_locations(nrs);

  // Loop over all broadening species and see if we can find them in abs_species.
  for (Index i = 0; i < nrs; ++i) {
    // Find associated internal species index (we do the lookup by index, not by name).
    const Species::Species isi = Species::fromShortName(ref_spec_names[i]);

    // Find position of broadening species isi in abs_species. The called
    // function returns -1 if not found, which is already the correct
    // treatment for this case that we also want here.
    ref_spec_locations[i] = find_first_species(abs_species, isi);
  }

  // The actual calculation
  // N_tot = sum (Nref_i *     p_i/p_0 * T0/T)
  //       = sum (Nref_i * vmr_i*p/p_0 * T0/T)
  //       = p/p_0 * T0/T *  sum (  Nref_i  * vmr_i)

  const Numeric ratioT = T0 / rtp_temperature;
  const Numeric ratiop = rtp_pressure / p0;

  Numeric ref_spec_vmr_sum = 0.;
  Numeric n = 0.;

  // Add up refractive species, where available:
  for (Index i = 0; i < nrs; ++i) {
    if (ref_spec_locations[i] >= 0) {
      // Add to VMR sum:
      ref_spec_vmr_sum += rtp_vmr[ref_spec_locations[i]];

      // refraction contribution (excluding the constant factor p/p_0 * T0/T):
      n += ref_n[i] * rtp_vmr[ref_spec_locations[i]];
    }
  }

  /*
  if ( abs(ref_spec_vmr_sum-1) > 0.1 )
      {
        std::ostringstream os;
        os << "Error: The total VMR of all your defined refractive\n"
             << "species is " << ref_spec_vmr_sum
             << ", more than 10% " << "different from 1.\n";
        throw std::runtime_error(os.str());
      }
  */

  // normalize refractive index with the considered total VMR:
  if (ref_spec_vmr_sum != 0) n /= ref_spec_vmr_sum;

  // now applying the constant factor p/p_0 * T0/T:
  n *= (ratioT * ratiop);

  refr_index_air += n;
  refr_index_air_group += n;
}

void complex_refr_indexWaterVisibleNIRHarvey98(GriddedField3& complex_refr_index,
                                const Vector& data_f_grid,
                                const Vector& data_t_grid,
                                const Vector& density_water,    //Gin
                                const Index& only_valid_range) {
  const Index N_f = data_f_grid.size();
  const Index N_t = data_t_grid.size();
  const Index N_d = density_water.size();

  if (N_d > 1)
    ARTS_USER_ERROR_IF(N_d != N_f, R"--(
density_water must be a Vector of size 1 or must be of the size of *data_t_grid*
)--")

  complex_refr_index.resize(N_f, N_t, 2);
  complex_refr_index.set_grid_name(0, "Frequency");
  complex_refr_index.set_grid(0, data_f_grid);
  complex_refr_index.set_grid_name(1, "Temperature");
  complex_refr_index.set_grid(1, data_t_grid);
  complex_refr_index.set_grid_name(2, "Complex");
  complex_refr_index.set_grid(2, ArrayOfString{"real", "imaginary"});

  Numeric refractive_index;
  Numeric density = density_water[0];

  for (Index i_f = 0; i_f < N_f; i_f++) {
    for (Index i_t = 0; i_t < N_t; i_t++) {
      if (density_water.size() > 1 && i_t > 0) density = density_water[i_t];

      refractive_index_water_and_steam_VisNIR(refractive_index,
                                              only_valid_range,
                                              data_f_grid[i_f],
                                              data_t_grid[i_t],
                                              density);

      complex_refr_index.data(i_f, i_t, 0) = refractive_index;
    }
  }
  complex_refr_index.data(joker, joker, 1) = 0.;
}

/*===========================================================================
  === WSMs for complex_refr_index
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void complex_refr_indexConstant(GriddedField3& complex_refr_index,
                                const Numeric& refr_index_real,
                                const Numeric& refr_index_imag) {
  complex_refr_index.resize(1, 1, 2);
  complex_refr_index.set_grid_name(0, "Frequency");
  complex_refr_index.set_grid(0, Vector(1, 0));
  complex_refr_index.set_grid_name(1, "Temperature");
  complex_refr_index.set_grid(1, Vector(1, 0));
  complex_refr_index.set_grid_name(2, "Complex");
  complex_refr_index.set_grid(2, ArrayOfString{"real", "imaginary"});

  complex_refr_index.data(joker, joker, 0) = refr_index_real;
  complex_refr_index.data(joker, joker, 1) = refr_index_imag;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void complex_refr_indexWaterLiebe93(GriddedField3& complex_refr_index,
                                    const Vector& f_grid,
                                    const Vector &t_grid) {
  if (min(t_grid) < 250) {
  }

  const Index nf = f_grid.size();
  const Index nt = t_grid.size();

  complex_refr_index.resize(nf, nt, 2);
  complex_refr_index.set_grid_name(0, "Frequency");
  complex_refr_index.set_grid(0, f_grid);
  complex_refr_index.set_grid_name(1, "Temperature");
  complex_refr_index.set_grid(1, t_grid);
  complex_refr_index.set_grid_name(2, "Complex");
  complex_refr_index.set_grid(2, ArrayOfString{"real", "imaginary"});

  Matrix complex_n;
  for (Index t = 0; t < nt; ++t) {
    complex_n_water_liebe93(complex_n, f_grid, t_grid[t]);
    complex_refr_index.data(joker, t, joker) = complex_n;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void complex_refr_indexIceMatzler06(GriddedField3& complex_refr_index,
                                    const Vector& f_grid,
                                    const Vector& t_grid) {
  const Index nf = f_grid.size();
  const Index nt = t_grid.size();

  // Frequency must be between 10MHz and 3THz
  const Numeric f_min = 10e6;
  const Numeric f_max = 3e12;
  chk_if_in_range(
      "min of complex_refr_index f_grid", min(f_grid), f_min, f_max);
  chk_if_in_range(
      "max of complex_refr_index f_grid", max(f_grid), f_min, f_max);

  // Temperature must be between 213.16 to 272.16 K
  const Numeric t_min = 20.;
  const Numeric t_max = 280.;
  chk_if_in_range(
      "min of complex_refr_index t_grid", min(t_grid), t_min, t_max);
  chk_if_in_range(
      "max of complex_refr_index t_grid", max(t_grid), t_min, t_max);

  complex_refr_index.resize(nf, nt, 2);
  complex_refr_index.set_grid_name(0, "Frequency");
  complex_refr_index.set_grid(0, f_grid);
  complex_refr_index.set_grid_name(1, "Temperature");
  complex_refr_index.set_grid(1, t_grid);
  complex_refr_index.set_grid_name(2, "Complex");
  complex_refr_index.set_grid(2, ArrayOfString{"real", "imaginary"});

  Matrix complex_n;
  for (Index i_t = 0; i_t < nt; ++i_t) {
    complex_n_ice_matzler06(complex_n, f_grid, t_grid[i_t]);
    complex_refr_index.data(joker, i_t, joker) = complex_n;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void complex_refr_indexTemperatureConstant(GriddedField3& complex_refr_index,
                                 const Vector& f_grid,
                                 const Vector& refr_index_real,
                                 const Vector& refr_index_imag,
                                 const Numeric& temperature) {
  chk_vector_length("f_grid","refr_index_real",f_grid,refr_index_real);
  chk_vector_length("f_grid","refr_index_imag",f_grid,refr_index_imag);

  const Index nf = f_grid.size();

  complex_refr_index.resize(nf, 1, 2);
  complex_refr_index.set_grid_name(0, "Frequency");
  complex_refr_index.set_grid(0, f_grid);
  complex_refr_index.set_grid_name(1, "Temperature");
  complex_refr_index.set_grid(1, Vector(1, temperature));
  complex_refr_index.set_grid_name(2, "Complex");
  complex_refr_index.set_grid(2, ArrayOfString{"real", "imaginary"});

  complex_refr_index.data(joker, 0, 0) = refr_index_real;
  complex_refr_index.data(joker, 0, 1) = refr_index_imag;
}

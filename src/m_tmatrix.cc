/*!
  \file   m_tmatrix.cc
  \author Oliver Lemke
  \date   2013-06-25

  \brief  T-Matrix related workspace methods.
*/

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include "arts_constants.h"
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "refraction.h"
#include "special_interp.h"
#include "tmatrix.h"

inline constexpr Numeric PI=Constant::pi;

/* Workspace method: Doxygen documentation will be auto-generated */
void diameter_maxFromDiameter_volume_equ(Numeric& diameter_max,
                                         Numeric& diameter_aspect_area_max,
                                         const String& shape,
                                         const Numeric& diameter_volume_equ,
                                         const Numeric& aspect_ratio) {
  const Numeric volume = (PI * pow(diameter_volume_equ, 3)) / 6.;

  if (shape == "spheroidal") {
    if (aspect_ratio < 1)  // prolate spheroid
    {
      //non-rotational axis (perpendicular to a)
      const Numeric b =
          pow((3. * volume) / (4. * PI * pow(aspect_ratio, 2)), 1. / 3.);
      diameter_max = 2. * b;
      diameter_aspect_area_max = 2. * b;
    } else  // oblate spheriod
    {
      //rotational axis
      const Numeric a = pow((3. * volume * aspect_ratio) / (4 * PI), 1. / 3.);
      diameter_max = 2. * a;
      diameter_aspect_area_max = 2. * a;
    }
  }

  else if (shape == "cylindrical") {
    //aspect_ratio=D/L
    const Numeric D = pow((volume * 4 * aspect_ratio) / PI, 1. / 3.);
    const Numeric L = D / aspect_ratio;
    diameter_max = pow(pow(D, 2) + pow(L, 2), 1. / 2.);
    diameter_aspect_area_max = std::max(D, std::pow(4 / PI * D * L, 1. / 2.));
  }

  else {
    std::ostringstream os;
    os << "Unknown particle shape: " << shape << "\n"
       << "Must be spheroidal or cylindrical";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void diameter_volume_equFromDiameter_max(Numeric& diameter_volume_equ,
                                         Numeric& volume,
                                         const String& shape,
                                         const Numeric& diameter_max,
                                         const Numeric& aspect_ratio) {
  if (shape == "spheroidal") {
    if (aspect_ratio < 1)  // prolate spheroid
    {
      const Numeric b = diameter_max / 2.;
      volume = (pow(b, 3.) * 4. * PI * pow(aspect_ratio, 2.)) / 3.;
    } else  // oblate spheriod
    {
      const Numeric a = diameter_max / 2.;
      volume = (pow(a, 3.) * 4. * PI) / (3. * aspect_ratio);
    }
  }

  else if (shape == "cylindrical") {
    volume = pow(diameter_max / pow((pow(aspect_ratio, 2.) + 1.), 1. / 2.), 3) *
             pow(aspect_ratio, 2.) * PI / 4.;
  }

  else {
    std::ostringstream os;
    os << "Unknown particle shape: " << shape << "\n"
       << "Must be spheroidal or cylindrical";
    throw std::runtime_error(os.str());
  }

  diameter_volume_equ = pow((6. * volume) / PI, 1. / 3.);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_singleTmatrix(SingleScatteringData& scat_data_single,
                             ScatteringMetaData& scat_meta_single,
                             const GriddedField3& complex_refr_index,
                             const String& shape,
                             const Numeric& diameter_volume_equ,
                             const Numeric& aspect_ratio,
                             const Numeric& mass,
                             const String& ptype,
                             const Vector& data_f_grid,
                             const Vector& data_t_grid,
                             const Vector& data_za_grid,
                             const Vector& data_aa_grid,
                             const Numeric& precision,
                             const String& cri_source,
                             const Index& ndgs,
                             const Index& robust,
                             const Index& quiet) {
  // Get internal coding for ptype
  scat_data_single.ptype = PTypeFromString(ptype);

  // Set description
  {
    std::ostringstream os;
    os << "T-matrix calculation for a " << shape << " particle, with "
       << "diameter_volume_equ = " << 1e6 * diameter_volume_equ << "um and "
       << "aspect ratio = " << aspect_ratio << ".";
    scat_data_single.description = os.str();
  }

  // Add grids to scat_data_single
  //
  scat_data_single.f_grid = data_f_grid;
  scat_data_single.T_grid = data_t_grid;

  if (scat_data_single.ptype == PTYPE_TOTAL_RND) {
    // tmatrix random orient requires equidistant angular grid. checking here
    // that given data_za_grid fulfills this requirement
    if (!(is_same_within_epsilon(data_za_grid[0], 0., 2 * DBL_EPSILON) &&
          is_same_within_epsilon(last(data_za_grid), 180., 2 * DBL_EPSILON))) {
      std::ostringstream os;
      os << "Zenith angle (=scattering angle) grid needs to include\n"
         << "0 deg and 180 deg as first and last grid points, respectively.\n"
         << "At least one of them does not fit.";
      throw std::runtime_error(os.str());
    }
    Index nza = data_za_grid.nelem();
    Numeric dza = 180. / ((Numeric)nza - 1.);
    for (Index iza = 1; iza < nza; iza++) {
      if (!(is_same_within_epsilon(
              data_za_grid[iza], (Numeric)iza * dza, 2 * DBL_EPSILON))) {
        std::ostringstream os;
        os << "Input zenith angle grid *data_za_grid* is required to be\n"
           << "equidistant for randomly oriented particles, but it is not.";
        throw std::runtime_error(os.str());
      }
    }
  }
  scat_data_single.za_grid = data_za_grid;

  if (scat_data_single.ptype == PTYPE_TOTAL_RND) {
    // in case of random orientation, azimuth grid should be empty. We just
    // set that here, ignoring whatever is in data_aa_grid.
    Vector empty_grid(0);
    scat_data_single.aa_grid = empty_grid;
  } else {
    // For azimuthally-random oriented particles, the azimuth angle grid must cover
    // 0-180 degrees.
    if (scat_data_single.ptype == PTYPE_AZIMUTH_RND &&
        data_aa_grid.nelem() == 0) {
      std::ostringstream os;
      os << "For ptype = \"azimuthally_random\""
         << " the azimuth angle grid can not be empty.";
      throw std::runtime_error(os.str());
    }
    if (scat_data_single.ptype == PTYPE_AZIMUTH_RND && data_aa_grid[0] != 0.) {
      std::ostringstream os;
      os << "For ptype = \"azimuthally_random\""
         << " the first value of the aa grid must be 0.";
      throw std::runtime_error(os.str());
    }

    if (scat_data_single.ptype == PTYPE_AZIMUTH_RND &&
        last(data_aa_grid) != 180.) {
      std::ostringstream os;
      os << "For ptype = \"azimuthally_random\""
         << " the last value of the aa grid must be 180.";
      throw std::runtime_error(os.str());
    }

    scat_data_single.aa_grid = data_aa_grid;
  }

  // Index coding for shape
  Index np;
  Numeric ar = aspect_ratio;
  if (shape == "spheroidal") {
    np = -1;
    if (aspect_ratio == 1)
      /*      do not throw error, but slightly increase aspect ratio to value
 *      recommended by original tmatrix code such that numerical issues are
 *      avoided.
        throw std::runtime_error( "For spheroidal particles, the aspect ratio "
                             "is not allowed to be exactly 1 (due to "
                             "numerical problems)." );
*/
      ar += 1e-6;
  } else if (shape == "cylindrical") {
    np = -2;
  } else {
    std::ostringstream os;
    os << "Unknown particle shape: " << shape << "\n"
       << "Must be \"spheroidal\" or \"cylindrical\".";
    throw std::runtime_error(os.str());
  }

  // Interpolate refractive index to relevant grids
  //
  const Index nf = data_f_grid.nelem();
  const Index nt = data_t_grid.nelem();
  //
  Tensor3 ncomp(nf, nt, 2);
  complex_n_interp(ncomp(joker, joker, 0),
                   ncomp(joker, joker, 1),
                   complex_refr_index,
                   "complex_refr_index",
                   data_f_grid,
                   data_t_grid);

  // Run T-matrix and we are ready (T-matrix takes size as volume equiv radius(!) )
  calcSingleScatteringDataProperties(scat_data_single,
                                     ncomp(joker, joker, 0),
                                     ncomp(joker, joker, 1),
                                     0.5 * diameter_volume_equ,
                                     np,
                                     ar,
                                     precision,
                                     ndgs,
                                     robust,
                                     quiet);

  // Meta data
  scat_meta_single.description =
      "Meta data for associated file with single scattering data.";
  scat_meta_single.source =
      "ARTS interface to T-matrix code by Mishchenko et al.";
  scat_meta_single.refr_index = cri_source;
  //
  Numeric diameter_max, area_max;
  diameter_maxFromDiameter_volume_equ(diameter_max,
                                      area_max,
                                      shape,
                                      diameter_volume_equ,
                                      aspect_ratio);
  //
  scat_meta_single.mass = mass;
  scat_meta_single.diameter_max = diameter_max;
  scat_meta_single.diameter_volume_equ = diameter_volume_equ;
  scat_meta_single.diameter_area_equ_aerodynamical = area_max;
}

void TMatrixTest() {
  tmatrix_tmd_test();
  tmatrix_ampld_test();
  calc_ssp_random_test();
  calc_ssp_fixed_test();
}

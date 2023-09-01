/**
 * @file   m_montecarlo.cc
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2003-06-19
 *
 * @brief  Workspace functions for the solution of cloud-box radiative transfer
 * by Monte Carlo methods.  All of these functions refer to 3D calculations
 *
 */
/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <ctime>
#include <fstream>
#include <stdexcept>
#include "arts_constants.h"
#include "arts_conversions.h"
#include <workspace.h>
#include "check_input.h"
#include "debug.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "mc_interp.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "refraction.h"
#include "rng.h"
#include "rte.h"
#include "special_interp.h"
#include "surf.h"
#include "xml_io.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric BOLTZMAN_CONST=Constant::boltzmann_constant;
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetGaussian(MCAntenna& mc_antenna,
                           //keyword arguments
                           const Numeric& za_sigma,
                           const Numeric& aa_sigma) {
  mc_antenna.set_gaussian(za_sigma, aa_sigma);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetGaussianByFWHM(MCAntenna& mc_antenna,
                                 //keyword arguments
                                 const Numeric& za_fwhm,
                                 const Numeric& aa_fwhm) {
  mc_antenna.set_gaussian_fwhm(za_fwhm, aa_fwhm);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetPencilBeam(MCAntenna& mc_antenna) {
  mc_antenna.set_pencil_beam();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MCSetSeedFromTime(Index& mc_seed) {
  mc_seed = (Index)time(NULL);
}

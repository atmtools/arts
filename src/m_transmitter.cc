/*===========================================================================
  ===  File description
  ===========================================================================*/

/**
  @file   m_transmitter.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2012-10-31

  @brief  Workspace functions related to transmitters and radiative transfer
  for transmitted signals.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts_constants.h"
#include "arts_conversions.h"
#include <workspace.h>
#include "geodetic.h"
#include "jacobian.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_complex.h"
#include "rte.h"
#include "rtepack.h"
#include "rtepack_multitype.h"
#include "sensor.h"
#include <cmath>
#include <stdexcept>

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterMultiplePol(Matrix& iy_transmitter,
                               const Vector& f_grid,
                               const ArrayOfIndex& instrument_pol) {
  const Index nf = f_grid.nelem();

  if (instrument_pol.nelem() != nf)
    throw std::runtime_error(
        "The length of *f_grid* and the number of elements "
        "in *instrument_pol* must be equal.");

  iy_transmitter.resize(nf, 4);

  for (Index i = 0; i < nf; i++) {
    stokes2pol(iy_transmitter(i, joker), instrument_pol[i], 1);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterSinglePol(Matrix& iy_transmitter,
                             const Vector& f_grid,
                             const ArrayOfIndex& instrument_pol) {
  const Index nf = f_grid.nelem();

  if (instrument_pol.nelem() != 1)
    throw std::runtime_error(
        "The number of elements in *instrument_pol* must be 1.");

  iy_transmitter.resize(nf, 4);

  stokes2pol(iy_transmitter(0, joker), instrument_pol[0], 1);

  for (Index i = 1; i < nf; i++) {
    iy_transmitter(i, joker) = iy_transmitter(0, joker);
  }
}

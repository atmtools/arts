/**
 * @file   m_physics.cc
 * @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
 * @date   2002-08-20
 *
 * @brief  Workspace methods of physical character.
 *
 * This file includes workspace methods for operations that have some
 * connection to basic physics. Example of methods are:  <br>
 * 1. Setting WSV to hold blackbody radiation. <br>
 * 2. Conversion to brightness temperature.
 *
 * These functions are listed in the doxygen documentation as entries of the
 * file auto_md.h.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts_constants.h"
#include <workspace.h>
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "physics_funcs.h"

inline constexpr Numeric TEMP_0_C=Constant::temperature_at_0c;
inline constexpr Numeric COSMIC_BG_TEMP=Constant::cosmic_microwave_background_temperature;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixPlanck(  // WS Output:
    Matrix& m,
    // WS Input:
    // WS Generic Input:
    const Vector& f,
    const Numeric& t) {
  const Index n = f.nelem();

  if (n == 0) throw std::runtime_error("The given frequency vector is empty.");

  m.resize(n, 4);
  m = 0;

  planck(m(joker, 0), f, t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixCBR(  // WS Output:
    Matrix& m,
    // WS Input:
    // WS Generic Input:
    const Vector& f) {
  MatrixPlanck(m, f, COSMIC_BG_TEMP);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixUnitIntensity(  // WS Output:
    Matrix& m,
    // WS Input:
    // WS Generic Input:
    const Vector& f) {
  const Index n = f.nelem();

  if (n == 0) throw std::runtime_error("The given frequency vector is empty.");

  m.resize(n, 4);
  m = 0;

  for (Index i = 0; i < n; i++) {
    m(i, 0) = 1.0;
  }
}

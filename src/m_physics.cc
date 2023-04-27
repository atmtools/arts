/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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

#include "arts.h"
#include "arts_constants.h"
#include "auto_md.h"
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"

inline constexpr Numeric TEMP_0_C=Constant::temperature_at_0c;
inline constexpr Numeric COSMIC_BG_TEMP=Constant::cosmic_microwave_background_temperature;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixCBR(  // WS Output:
    Matrix& m,
    // WS Input:
    const Index& stokes_dim,
    // WS Generic Input:
    const Vector& f,
    const Verbosity&) {
  const Index n = f.nelem();

  if (n == 0) throw runtime_error("The given frequency vector is empty.");

  m.resize(n, stokes_dim);
  m = 0;

  planck(m(joker, 0), f, COSMIC_BG_TEMP);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixPlanck(  // WS Output:
    Matrix& m,
    // WS Input:
    const Index& stokes_dim,
    // WS Generic Input:
    const Vector& f,
    const Numeric& t,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  const Index n = f.nelem();

  if (n == 0) throw runtime_error("The given frequency vector is empty.");

  out2 << "  Setting blackbody radiation for a temperature of " << t << " K.\n";

  m.resize(n, stokes_dim);
  m = 0;

  planck(m(joker, 0), f, t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixUnitIntensity(  // WS Output:
    Matrix& m,
    // WS Input:
    const Index& stokes_dim,
    // WS Generic Input:
    const Vector& f,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  const Index n = f.nelem();

  if (n == 0) throw runtime_error("The given frequency vector is empty.");

  out2 << "  Setting unpolarised radiation with an intensity of 1.\n";

  m.resize(n, stokes_dim);
  m = 0;

  for (Index i = 0; i < n; i++) {
    m(i, 0) = 1.0;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void water_p_eq_fieldMK05(Tensor3& water_p_eq_field,
                          const AtmField& atm_field,
                          const Index& only_liquid,
                          const Verbosity&) {
ARTS_USER_ERROR_IF(not atm_field.regularized, "Must have regular grid atmospheric field")
const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();

  const Index n1 = t_field.npages();
  const Index n2 = t_field.nrows();
  const Index n3 = t_field.ncols();

  water_p_eq_field.resize(n1, n2, n3);

  for (Index i = 0; i < n1; i++) {
    for (Index j = 0; j < n2; j++) {
      for (Index k = 0; k < n3; k++) {
        const Numeric t = t_field(i, j, k);
        if (t > TEMP_0_C || only_liquid) {
          water_p_eq_field(i, j, k) =
              exp(54.842763 - 6763.22 / t - 4.21 * log(t) + 0.000367 * t +
                  tanh(0.0415 * (t - 218.8)) *
                      (53.878 - 1331.22 / t - 9.44523 * log(t) + 0.014025 * t));
        } else {
          water_p_eq_field(i, j, k) =
              exp(9.550426 - 5723.265 / t + 3.53068 * log(t) - 0.00728332 * t);
        }
      }
    }
  }
}

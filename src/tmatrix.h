/* Copyright (C) 2013 Oliver Lemke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*!
  \file   tmatrix.h
  \author Oliver Lemke
  \date   2013-06-25

  \brief  Declarations for the T-Matrix interface.
*/

#ifndef tmatrix_h
#define tmatrix_h

#include "messages.h"


/** T-Matrix validation test.

 Executes the standard test included with the double precision T-Matrix code
 for nonspherical particles in a fixed orientation.
 Should give the same as running the 3rdparty/tmatrix/tmatrix_ampld executable.

 \author Oliver Lemke
 */
void tmatrix_ampld_test(const Verbosity& verbosity);


/** T-Matrix validation test.

 Executes the standard test included with the double precision T-Matrix code
 for randomly oriented nonspherical particles.
 Should give the same as running the 3rdparty/tmatrix/tmatrix_tmd executable.

 \author Oliver Lemke
 */
void tmatrix_tmd_test(const Verbosity& verbosity);


/** Single scattering properties calculation for randomly oriented particles.

 Two cases are calculated. One with oblate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types

params = {'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'aspect_ratio': 1.000001}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 And one with prolate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types

params = {'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'aspect_ratio': 0.7}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 \author Oliver Lemke
 */
void calc_ssp_random_test(const Verbosity& verbosity);

/** Single scattering properties calculation for particles with fixed orientation.

 Two cases are calculated. One with oblate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types
from PyARTS import constants

params = {'ptype': constants.PARTICLE_TYPE_HORIZ_AL,
          'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'aspect_ratio': 1.000001}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 And one with prolate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types
from PyARTS import constants

params = {'ptype': constants.PARTICLE_TYPE_HORIZ_AL,
          'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'aspect_ratio': 0.7}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 \author Oliver Lemke
 */
void calc_ssp_fixed_test(const Verbosity& verbosity);

#endif //  tmatrix_h


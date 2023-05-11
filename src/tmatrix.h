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

#include "optproperties.h"

/** Calculate SingleScatteringData properties.

 Port of calc_SSP function from PyARTS.

 \param[in,out] ssd         As input all grids must be set (f, T, za and aa).
                            These grids are given by ssd are used to calculate
                            pha_mat_data, ext_mat_data and abs_vec_data.
                            Also ssd.ptype must be set as input.
                            The output ssd has remaining fields set.
 \param[in] ref_index_real  Vector with real parts of refractive index
                            Number of rows must match elements in ssd.f_grid
                            Number of cols must match elements in ssd.T_grid
 \param[in] ref_index_imag  Vector with imaginary parts of refractive index
 \param[in] equiv_radius    equivalent volume radius [micrometer]
 \param[in] np              Particle type (-1 for spheroid, -2 for cylinder)
 \param[in] axial_ratio     Axial ratio of particles
 \param[in] precision       Accuracy of the computations
 \param[in] ndgs            The number of division points in computing
                            integrals over the particle surface.

 \author Oliver Lemke
 */
void calcSingleScatteringDataProperties(SingleScatteringData& ssd,
                                        ConstMatrixView ref_index_real,
                                        ConstMatrixView ref_index_imag,
                                        const Numeric equiv_radius = 200,
                                        const Index np = -1,
                                        const Numeric axial_ratio = 1.000001,
                                        const Numeric precision = 0.001,
                                        const Index ndgs = 2,
                                        const Index robust = 0,
                                        const Index quiet = 1);

/** T-Matrix validation test.

 Executes the standard test included with the double precision T-Matrix code
 for nonspherical particles in a fixed orientation.
 Should give the same as running the 3rdparty/tmatrix/tmatrix_ampld executable.

 \author Oliver Lemke
 */
void tmatrix_ampld_test();

/** T-Matrix validation test.

 Executes the standard test included with the double precision T-Matrix code
 for randomly oriented nonspherical particles.
 Should give the same as running the 3rdparty/tmatrix/tmatrix_tmd executable.

 \author Oliver Lemke
 */
void tmatrix_tmd_test();

/** Single scattering properties calculation for randomly oriented particles.

 Two cases are calculated. One with oblate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types

params = {'ptype': constants.PTYPE_AZIMUTH_RND,
          'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'mrr': numpy.array([[1.78031135, 1.78150475], [1.78037238, 1.78147686]]),
          'mri': numpy.array([[0.00278706, 0.00507565], [0.00287245, 0.00523012]]),
          'axial_ratio': 1.5}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 And one with prolate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types

params = {'ptype': constants.PTYPE_AZIMUTH_RND,
          'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'mrr': numpy.array([[1.78031135, 1.78150475], [1.78037238, 1.78147686]]),
          'mri': numpy.array([[0.00278706, 0.00507565], [0.00287245, 0.00523012]]),
          'axial_ratio': 0.7}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 \author Oliver Lemke
 */
void calc_ssp_random_test();

/** Single scattering properties calculation for particles with fixed orientation.

 Two cases are calculated. One with oblate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types
from PyARTS import constants

params = {'ptype': constants.PTYPE_AZIMUTH_RND,
          'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'mrr': numpy.array([[1.78031135, 1.78150475], [1.78037238, 1.78147686]]),
          'mri': numpy.array([[0.00278706, 0.00507565], [0.00287245, 0.00523012]]),
          'axial_ratio': 1.5}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 And one with prolate particles which is equivalent to
 the following PyARTS case:
 <pre>
from PyARTS import arts_types
from PyARTS import constants

params = {'ptype': constants.PTYPE_AZIMUTH_RND,
          'f_grid': [230e9, 240e9],
          'T_grid': [220, 250],
          'za_grid': numpy.arange(0, 181, 10),
          'aa_grid': numpy.arange(0, 181, 10),
          'equiv_radius': 200, # equivalent volume radius
          'NP':-1, # -1 for spheroid, -2 for cylinder, positive for chebyshev
          'phase':'ice',
          'mrr': numpy.array([[1.78031135, 1.78150475], [1.78037238, 1.78147686]]),
          'mri': numpy.array([[0.00278706, 0.00507565], [0.00287245, 0.00523012]]),
          'axial_ratio': 0.7}
s = arts_types.SingleScatteringData(params)
s.calc()
 </pre>

 \author Oliver Lemke
 */
void calc_ssp_fixed_test();

#endif  //  tmatrix_h

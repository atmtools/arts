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
  \file   m_tmatrix.cc
  \author Oliver Lemke
  \date   2013-06-25

  \brief  T-Matrix related workspace methods.
*/

#include <stdexcept>
#include "messages.h"
#include "tmatrix.h"

void TMatrixTest(const Verbosity& verbosity)
{
    tmatrix_tmd_test(verbosity);
    tmatrix_ampld_test(verbosity);
    calc_ssp_random_test(verbosity);
    calc_ssp_fixed_test(verbosity);
}

void single_scattering_dataCalcTMatrixTest(// WS Output:
                                           SingleScatteringData& sdd,
                                           // WS Generic input:
                                           const String& p_type,
                                           const Vector& f_grid,
                                           const Vector& T_grid,
                                           const Vector& za_grid,
                                           const Vector& aa_grid,
                                           const Matrix& ref_index_real,
                                           const Matrix& ref_index_imag,
                                           const Numeric& equiv_radius,
                                           const Index& np,
                                           const String& phase _U_,
                                           const Numeric& aspect_ratio,
                                           const Numeric& precision,
                                           const Verbosity&)
{
    sdd.f_grid = f_grid;
    sdd.T_grid = T_grid;
    sdd.za_grid = za_grid;
    sdd.aa_grid = aa_grid;

    if (p_type == "MACROS_ISO")
        sdd.ptype = PARTICLE_TYPE_MACROS_ISO;
    else if (p_type == "HORIZ_AL")
        sdd.ptype = PARTICLE_TYPE_HORIZ_AL;
    else
    {
        ostringstream os;
        os << "Unknown particle type: " << p_type << "\n"
        << "Must be MACROS_ISO or HORIZ_AL";
        throw std::runtime_error(os.str());
    }

    calcSingleScatteringDataProperties(sdd,
                                       ref_index_real,
                                       ref_index_imag,
                                       equiv_radius,
                                       np,
                                       phase,
                                       aspect_ratio,
                                       precision);
}

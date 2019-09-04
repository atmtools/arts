/* Copyright (C) 2019
   Richard Larsson
                            
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
 * @file radiation_field.cc
 * @author Richard Larsson
 * @date 2019-09-04
 * 
 * @brief Radiation field calculations
 */

#include "radiation_field.h"
#include "sorting.h"

void error_in_integrate(const String& error_msg,
                        const Numeric& value_that_should_be_unity) {
  if (std::abs(value_that_should_be_unity - 1.0) > 1e-4) {
    std::ostringstream os;
    os << "Failure in normalization:\n" << error_msg << "\n";
    throw std::runtime_error(os.str());
  }
}

Numeric test_integrate_convolved(const Eigen::Ref<Eigen::VectorXcd> F,
                                 const Vector& f) {
  Numeric val = 0.0;

  const Index n = f.nelem();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (f[i + 1] - f[i]) * (F[i].real() + F[i + 1].real());

  return val;  // Should return 1.0 for normalized F
}

Numeric test_integrate_zenith(const Vector& cosza,
                              const Array<Index>& sorted_index) {
  Numeric val = 0.0;

  const Index n = cosza.nelem();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (cosza[sorted_index[i]] - cosza[sorted_index[i + 1]]);

  return val;  // Should return 1.0 for normalized cosza
}

Numeric integrate_convolved(const RadiationVector& I,
                            const Eigen::VectorXcd& F,
                            const Vector& f) {
  Numeric val = 0.0;

  const Index n = f.nelem();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (f[i + 1] - f[i]) *
           (I(i, 0) * F[i].real() + I(i + 1, 0) * F[i + 1].real());

  return val;
}

Numeric integrate_convolved(const TransmissionMatrix& T,
                            const Eigen::VectorXcd& F,
                            const Vector& f) {
  Numeric val = 0.0;

  const Index n = f.nelem();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (f[i + 1] - f[i]) *
           (T(i, 0, 0) * F[i].real() + T(i + 1, 0, 0) * F[i + 1].real());

  return 1.0 - val;
}

Numeric integrate_zenith(const VectorView j,
                         const Vector& cosza,
                         const Array<Index>& sorted_index) {
  Numeric val = 0.0;

  const Index n = cosza.nelem();
  for (Index i = 0; i < n - 1; i++)
    val += 0.25 * (cosza[sorted_index[i]] - cosza[sorted_index[i + 1]]) *
           (j[sorted_index[i]] + j[sorted_index[i + 1]]);

  return val;
}

Index grid_index_from_gp(const GridPos& gp) {
  if (gp.fd[1] == 1.0)
    return gp.idx;
  else
    return gp.idx + 1;
}

void sorted_index_of_ppath_field(ArrayOfArrayOfIndex& sorted_index,
                                 ArrayOfVector& cos_zenith_angles,
                                 const ArrayOfPpath& ppath_field) {
  extern const Numeric DEG2RAD;

  Index nalt = 0;
  for (auto& path : ppath_field)
    for (auto& gp : path.gp_p) nalt = std::max(nalt, grid_index_from_gp(gp));
  nalt += 1;

  // Get the data, which is of unknown length
  Array<ArrayOfNumeric> zeniths_array(nalt, ArrayOfNumeric(0));
  for (auto& path : ppath_field) {
    for (Index ip = 0; ip < path.np; ip++) {
      const Index ind = grid_index_from_gp(path.gp_p[ip]);
      zeniths_array[ind].push_back(path.los(ip, 0));
    }
  }

  // Finalize sorting
  sorted_index.resize(nalt);
  cos_zenith_angles.resize(nalt);
  for (Index i = 0; i < nalt; i++) {
    Vector& data = cos_zenith_angles[i];
    data.resize(zeniths_array[i].nelem());

    for (Index j = 0; j < data.nelem(); j++) data[j] = zeniths_array[i][j];
    get_sorted_indexes(sorted_index[i], data);

    for (Index j = 0; j < data.nelem(); j++)
      data[j] = std::cos(DEG2RAD * data[j]);
  }
}

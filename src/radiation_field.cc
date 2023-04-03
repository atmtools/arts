/**
 * @file radiation_field.cc
 * @author Richard Larsson
 * @date 2019-09-04
 * 
 * @brief Radiation field calculations
 */

#include "radiation_field.h"
#include "arts_conversions.h"
#include "sorting.h"

void error_in_integrate(const String& error_msg,
                        const Numeric& value_that_should_be_unity) {
  ARTS_USER_ERROR_IF (std::abs(value_that_should_be_unity - 1.0) > 1e-4,
    "Failure in normalization:\n", error_msg, "\n")
}

Numeric test_integrate_convolved(const Eigen::Ref<Eigen::VectorXcd> &F,
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

Numeric integrate_zenith(const ConstVectorView& j,
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
  return gp.idx + 1;
}

void sorted_index_of_ppath_field(ArrayOfArrayOfIndex& sorted_index,
                                 ArrayOfVector& cosza,
                                 const ArrayOfPpath& ppath_field) {
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
  cosza.resize(nalt);
  for (Index i = 0; i < nalt; i++) {
    Vector& data = cosza[i];
    data.resize(zeniths_array[i].nelem());

    for (Index j = 0; j < data.nelem(); j++) data[j] = zeniths_array[i][j];
    get_sorted_indexes(sorted_index[i], data);

    for (Index j = 0; j < data.nelem(); j++)
      data[j] = Conversion::cosd(data[j]);
  }
}
